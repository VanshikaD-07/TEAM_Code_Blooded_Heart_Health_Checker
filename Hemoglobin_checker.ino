#include <Wire.h>
#include "MAX30105.h"
#include <Adafruit_GFX.h>
#include <Adafruit_SSD1306.h>
#include <math.h>

#define SDA_PIN    21
#define SCL_PIN    22
#define BUTTON_PIN 4

MAX30105 sensor;
Adafruit_SSD1306 oled(128, 64, &Wire, -1);

// ---------------- Sensor settings ----------------
const byte AMP_RED = 60;           // 0-255. Tune so red/IR DC are ~50k-200k
const byte AMP_IR  = 60;
const byte SAMPLE_AVG  = 4;
const byte LED_MODE    = 2;        // 2 = red + IR only
const int  SAMPLE_RATE = 400;      // 400 Hz / 4 averaging = 100 Hz per channel
const int  PULSE_WIDTH = 411;
const int  ADC_RANGE   = 16384;

const int   FS      = 100;
const float DT_MS   = 1000.0f / FS;
const int   WINDOW_SEC = 30;                 // keep <= 45 (RAM)
const int   N          = FS * WINDOW_SEC;
const int   SETTLE_SAMPLES = 2 * FS;         // ignore first 2 s after finger placed
const int   EDGE = FS;                       // trim 1 s at each end after filtering
const long  FINGER_THRESH = 50000;
const int   MAX_BEATS = 120;

// ---------------- Quality-gate thresholds ----------------
const float DC_MIN = 30000, DC_MAX = 250000;
const float PI_MIN_PCT   = 0.05f;   // IR perfusion index, %
const float SNR_MIN      = 3.0f;    // pulse RMS / high-frequency noise RMS
const float CORR_MIN     = 0.85f;   // red vs IR waveform similarity
const float DRIFT_MAX_PCT = 5.0f;   // baseline wander (pressure/motion)
const int   NN_MIN       = 15;      // minimum valid beat-to-beat intervals
const float KEEP_FRAC_MIN = 0.80f;  // fraction of intervals that are plausible
const int   WAVE_MIN     = 8;       // minimum beats usable for morphology
const int   BA_MIN       = 5;

enum { F_DC = 1, F_PERF = 2, F_SNR = 4, F_CORR = 8, F_DRIFT = 16,
       F_BEATS = 32, F_INTERVALS = 64, F_WAVE = 128 };

// ---------------- Hb model (EDIT AFTER CALIBRATION) ----------------
// Hb = baseline(sex) + sum_k WEIGHT[k] * (feature_k - FEAT_MEAN[k]) / FEAT_SD[k]
// Until you fit these from volunteer data, keep MODEL_CALIBRATED = false.
enum { SEX_UNSET = 0, SEX_MALE = 1, SEX_FEMALE = 2 };
const float BASE_HB_MALE   = 14.5f;   // replace with mean lab Hb of YOUR male volunteers
const float BASE_HB_FEMALE = 12.5f;   // replace with mean lab Hb of YOUR female volunteers
const float ANEMIA_CUTOFF_MALE   = 13.0f;   // WHO adult, g/dL
const float ANEMIA_CUTOFF_FEMALE = 12.0f;
const bool  MODEL_CALIBRATED = false;
const float MODEL_MAE = 0.0f;               // cross-validated error (g/dL) from your Python fit

const int K = 12;
const char *FEAT_NAMES[K] = { "dc_ratio", "R_red", "R_ir", "R_ratio", "pi_pct", "hr_bpm",
                              "rmssd_ms", "rise_ms", "width50_ms", "area_norm", "ba_ratio", "delay_ms" };
float FEAT_MEAN[K] = {0,0,0,0,0,0,0,0,0,0,0,0};
float FEAT_SD[K]   = {1,1,1,1,1,1,1,1,1,1,1,1};
float WEIGHT[K]    = {0,0,0,0,0,0,0,0,0,0,0,0};

const bool PRETTY_PRINT = true;
const bool LOG_FAILED_ROWS = false;

// ---------------- Buffers ----------------
uint32_t bufRed[N], bufIR[N];
float fRed[N], fIr[N], scratch[N];
float gPos[MAX_BEATS], gRaw[MAX_BEATS];
bool  gValid[MAX_BEATS];
float gRise[MAX_BEATS], gW50[MAX_BEATS], gArea[MAX_BEATS], gBA[MAX_BEATS];

// ---------------- State ----------------
int      count = 0, settle = 0;
bool     needRemoval = false, waitShown = false;
uint32_t holdUntil = 0, lastButton = 0;
int      subjectSex = SEX_UNSET;
String   subjectId = "none";
float    refHb = NAN;
String   cmdBuf;

struct Result {
  uint16_t failMask;
  bool  ok;
  float redDc, irDc, dcRatio, redAc, irAc, Rred, Rir, Rratio, spo2, pi;
  float hr, sdnn, rmssd;
  int   nBeats, nNN;
  float rise, width50, areaNorm, ba, delayMs;
  float snr, corr, drift;
  float prior, hb;
};

// =====================================================================
// Function declared
// =====================================================================
float medianOf(const float *v, int n) {
  if (n <= 0) return NAN;
  float t[MAX_BEATS];
  if (n > MAX_BEATS) n = MAX_BEATS;
  for (int i = 0; i < n; i++) t[i] = v[i];
  for (int i = 1; i < n; i++) {
    float key = t[i]; int j = i - 1;
    while (j >= 0 && t[j] > key) { t[j + 1] = t[j]; j--; }
    t[j + 1] = key;
  }
  return (n & 1) ? t[n / 2] : 0.5f * (t[n / 2 - 1] + t[n / 2]);
}

char sexChar() { return subjectSex == SEX_MALE ? 'M' : (subjectSex == SEX_FEMALE ? 'F' : '?'); }

// Detrend (1 s centred moving average), invert (more blood = less light, so
// systolic peaks point up), smooth (5-pt). Returns pulse RMS (the "AC" part).
float prepareSignal(const uint32_t *x, float *f, float dc, float &driftPct, float &hfRms) {
  const int HALF = FS / 2;
  float bMin = 1e12f, bMax = -1e12f;
  for (int i = 0; i < N; i++) {
    int lo = i - HALF < 0 ? 0 : i - HALF;
    int hi = i + HALF > N - 1 ? N - 1 : i + HALF;
    double s = 0;
    for (int k = lo; k <= hi; k++) s += x[k];
    float base = (float)(s / (hi - lo + 1));
    scratch[i] = -((float)x[i] - base);
    if (i >= EDGE && i < N - EDGE) {
      if (base < bMin) bMin = base;
      if (base > bMax) bMax = base;
    }
  }
  for (int i = 0; i < N; i++) {
    int lo = i - 2 < 0 ? 0 : i - 2;
    int hi = i + 2 > N - 1 ? N - 1 : i + 2;
    float s = 0;
    for (int k = lo; k <= hi; k++) s += scratch[k];
    f[i] = s / (hi - lo + 1);
  }
  double sq = 0, hf = 0;
  for (int i = EDGE; i < N - EDGE; i++) {
    sq += (double)f[i] * f[i];
    float m = 0;
    for (int k = -2; k <= 2; k++) m += (float)x[i + k];
    m /= 5.0f;
    float d = (float)x[i] - m;
    hf += (double)d * d;
  }
  int n = N - 2 * EDGE;
  hfRms = (float)sqrt(hf / n);
  driftPct = 100.0f * (bMax - bMin) / dc;
  return (float)sqrt(sq / n);
}
