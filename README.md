# 🌊 LoRaWAN-Based Marine Monitoring System (MATLAB)

## 🚀 Overview

This project implements a **complete end-to-end LoRa physical layer simulation** for long-range marine communication, enhanced with a **LoRaWAN Adaptive Data Rate (ADR) algorithm**.

It models a real-world **Arabian Sea buoy system**, transmitting environmental data over a **15 km wireless link** under dynamic ocean conditions.

---

## 🧠 Key Highlights

### 📡 Physical Layer (LoRa PHY)

* Chirp Spread Spectrum (CSS) modulation
* LoRa packet formation
* Hamming coding + Gray mapping
* BER vs SNR (Monte Carlo simulation)

---

### ⚡ Hardware Impairments

* Carrier Frequency Offset (CFO)
* Phase noise
* Clock drift

---

### 🌊 Channel Modeling

* AWGN
* Rayleigh fading
* Rician fading
* Log-normal shadowing
* Time-varying SNR (sea-state + monsoon effects)

---

## 🚀 Adaptive Data Rate (ADR) — Advanced Feature ⭐

This project implements a **full LoRaWAN ADR algorithm** with realistic network behavior:

### ✔ Features

* SNR history-based decision (last 20 packets)
* Dynamic Spreading Factor (SF7–SF12)
* Adaptive transmit power control
* Margin-based optimization

### ✔ Realistic Scenario

* 24-hour buoy deployment (288 packets)
* Time-varying SNR profile:

  * Morning: calm sea (high SNR)
  * Evening: monsoon dip (−8 dB)
  * Night: recovery phase

---

## 📊 Performance Comparison

| Metric          | Fixed SF10 | ADR      |
| --------------- | ---------- | -------- |
| Packet Delivery | ~X%        | ~X%      |
| Energy Usage    | High       | Reduced  |
| Airtime         | High       | Reduced  |
| Battery Life    | Lower      | Extended |

👉 ADR dynamically:

* Reduces SF → faster transmission (good channel)
* Increases SF → reliable communication (poor channel)
* Lowers power → saves battery

---

## ⚙️ System Architecture

Sensor → Encoding → LoRa Packet → Chirp Modulation → Channel → Receiver → Decoding → ADR Optimization

---

## 📁 Project Structure

```id="z8q1vp"
src/        → MATLAB source code (PHY + ADR)
results/    → Plots (BER, ADR performance, SNR profile)
docs/       → Report / diagrams
data/       → Simulation data
```

---

## ▶️ How to Run

Run scripts in order:

```id="ptj7q9"
Chunk1_Setup
Chunk2_Modulation
Chunk3_Channel
Chunk4_BER_Decode
chunk5_system_analysis
Chunk6_ADR
```

---

## 🛠 Tools & Technologies

* MATLAB
* Wireless Communication Systems
* Signal Processing
* LPWAN (LoRaWAN)

---

## 🌍 Applications

* Marine monitoring systems
* IoT-based long-range communication
* Smart ocean sensing
* Wireless PHY research

---

## 📈 Future Improvements

* SDR-based real hardware implementation
* LoRaWAN MAC layer integration
* Network-level multi-node simulation

---

## 👤 Author

**Unnimaya**

---

## ⭐ If you like this project

Give it a ⭐ on GitHub!

