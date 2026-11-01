"""Save the ESP32's CSV rows to a file on your laptop.

Setup (once):   pip install pyserial
Run:            python log_serial.py COM5            (Windows; check the port in Arduino IDE)
                python log_serial.py /dev/ttyUSB0    (Linux)
                python log_serial.py /dev/cu.usbserial-0001   (macOS)
Close the Arduino Serial Monitor first: only one program can use the port.

You can type the firmware commands here and press Enter:
    m / f   set sex      r13.2  set reference Hb      iA01  set subject ID      x  clear Hb
Quit with Ctrl+C. Rows are appended to hb_data.csv in the folder you run this from.
"""
import sys
import threading
import time
from pathlib import Path

import serial  # pip install pyserial

HEADER = ("id,sex,ref_hb,ts_ms,quality_ok,fail_mask,hr_bpm,rmssd_ms,sdnn_ms,spo2,pi_pct,"
          "red_dc,ir_dc,dc_ratio,R_red,R_ir,R_ratio,rise_ms,width50_ms,area_norm,ba_ratio,"
          "delay_ms,snr,corr,drift_pct,n_beats,n_nn,hb_est")
N_COMMAS = HEADER.count(",")
OUT = Path("hb_data.csv")


def main():
    if len(sys.argv) < 2:
        sys.exit("Usage: python log_serial.py <serial-port>   e.g. COM5")
    ser = serial.Serial(sys.argv[1], 115200, timeout=1)
    time.sleep(2)  # the ESP32 resets when the port opens

    if not OUT.exists() or OUT.stat().st_size == 0:
        OUT.write_text(HEADER + "\n")
    print(f"Logging data rows to: {OUT.resolve()}\n")

    def send_commands():
        for line in sys.stdin:
            ser.write((line.strip() + "\n").encode())

    threading.Thread(target=send_commands, daemon=True).start()

    saved = 0
    try:
        while True:
            raw = ser.readline().decode("utf-8", errors="replace").rstrip()
            if not raw:
                continue
            print(raw)
            if raw.count(",") == N_COMMAS and raw != HEADER:
                with OUT.open("a") as f:
                    f.write(raw + "\n")
                saved += 1
                print(f"   >>> row {saved} saved to {OUT.name}")
    except KeyboardInterrupt:
        print(f"\nStopped. {saved} rows saved in {OUT.resolve()}")


if __name__ == "__main__":
    main()
