# 🌐 Remote Access & Deployment Guide

This guide details methods for accessing a PC- or server-hosted **CLAIREscope** instance remotely from laptops, tablets, mobile devices, and collaborating research groups.

---

## 1. NordVPN Meshnet (Recommended for Secure Multi-Device / Lab Access)

**NordVPN Meshnet** creates an encrypted peer-to-peer virtual private network connecting your devices and collaborators without exposing open ports or requiring complex router configuration.

### Setup Instructions
1. **Enable Meshnet on Host Machine (Windows / Server)**:
   - Open the **NordVPN app** $\rightarrow$ Click the **Meshnet icon** in the left sidebar $\rightarrow$ Toggle **Meshnet ON**.
   - Note the assigned **Meshnet IP** (e.g. `100.64.x.x`) or **Nord Name** (e.g. `claire-desktop.nord`).

2. **Configure Windows WSL2 Port Forwarding** (Run once in PowerShell as Administrator):
   ```powershell
   netsh interface portproxy add v4tov4 listenport=8501 listenaddress=0.0.0.0 connectport=8501 connectaddress=127.0.0.1
   ```

3. **Access from Connected Devices**:
   - Turn on Meshnet in the NordVPN client on your client device (phone, iPad, laptop) or send a Meshnet invitation link to collaborators.
   - Access the web interface in any browser:
     ```text
     http://<your-meshnet-ip>:8501
     # or
     http://<your-nord-name>:8501
     ```

---

## 2. Cloudflare Tunnels (`cloudflared`) (Instant HTTPS for Collaborators)

Generates an ephemeral, encrypted public HTTPS link (`https://xxxx.trycloudflare.com`) with zero network or firewall modifications.

### Quick Start
```bash
# Install and run cloudflared ephemeral tunnel pointing to CLAIREscope port
cloudflared tunnel --url http://localhost:8501
```
Share the output URL with collaborators for instant browser access without client VPN installations.

---

## 3. Local Area Network (LAN / Wi-Fi) Access

For devices connected to the same local Wi-Fi or lab network:
```text
http://<host-local-ip>:8501
# e.g., http://192.168.146.133:8501
```

---

## 4. Permanent Cloud Hosting (Streamlit Community Cloud)

For public demonstration datasets:
1. Push dataset or demo scripts to GitHub (`https://github.com/ccneko/CLAIREscope`).
2. Log into [share.streamlit.io](https://share.streamlit.io) and deploy repository `ccneko/CLAIREscope` with main file `app.py`.
3. Provides a permanent 24/7 URL (e.g. `https://clairescope.streamlit.app`).
