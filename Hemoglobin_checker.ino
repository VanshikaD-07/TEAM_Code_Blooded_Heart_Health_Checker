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

