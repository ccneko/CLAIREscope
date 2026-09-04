# 🌐 Remote Access & Deployment Guide

This guide details methods for accessing a PC- or server-hosted **CLAIREscope** instance remotely from laptops, tablets, mobile devices, and collaborating research groups.

---

## 1. Local Area Network (LAN / Wi-Fi) Access

For devices connected to the same local Wi-Fi or laboratory network:

### Access Instructions
1. Find your host machine's local IPv4 address (e.g. \192.168.x.x\ via \ipconfig\ on Windows or \hostname -I\ in Linux/WSL).
2. Ensure the Streamlit server is running with \--server.address 0.0.0.0\:
   \\ash
   streamlit run app.py --server.port 8501 --server.address 0.0.0.0
   \3. Open any browser on your connected phone, tablet, or laptop and navigate to:
   \\	ext
   http://<host-local-ip>:8501
   # e.g., http://192.168.1.100:8501
   \
---

## 2. Cloudflare Tunnels (\cloudflared\) (Instant HTTPS for Collaborators)

Generates an ephemeral or persistent, encrypted public HTTPS link (\https://xxxx.trycloudflare.com\) with zero router, firewall, or port-forwarding modifications.

### Quick Start
\\ash
# Install and run cloudflared ephemeral tunnel pointing to CLAIREscope port
cloudflared tunnel --url http://localhost:8501
\Share the output URL with collaborators for instant secure browser access without client VPN installations.

---

## 3. NordVPN Meshnet (Recommended for Secure Multi-Device / Lab Access)

**NordVPN Meshnet** creates an encrypted peer-to-peer virtual private network connecting your personal devices and collaborating researchers without exposing open ports.

### Setup Instructions
1. **Enable Meshnet on Host Machine (Windows / Server)**:
   - Open the **NordVPN app** $\rightarrow$ Click the **Meshnet icon** in the left sidebar $\rightarrow$ Toggle **Meshnet ON**.
   - Note the assigned **Meshnet IP** (e.g. \100.64.x.x\) or **Nord Name** (e.g. \claire-desktop.nord\).

2. **Configure Windows WSL2 Port Forwarding** (Run once in PowerShell as Administrator if running inside WSL):
   \\powershell
   netsh interface portproxy add v4tov4 listenport=8501 listenaddress=0.0.0.0 connectport=8501 connectaddress=127.0.0.1
   \
3. **Access from Connected Devices**:
   - Turn on Meshnet in the NordVPN client on your client device (phone, iPad, laptop) or send a Meshnet invitation link to collaborators.
   - Access the web interface in any browser:
     \\	ext
     http://<your-meshnet-ip>:8501
     # or
     http://<your-nord-name>:8501
     \
---

## 4. Permanent Cloud Hosting (Streamlit Community Cloud)

For public demonstration datasets:
1. Push demo code or dataset configuration to GitHub (\https://github.com/ccneko/CLAIREscope\).
2. Log into [share.streamlit.io](https://share.streamlit.io) and deploy repository \ccneko/CLAIREscope\ with main file \pp.py\.
3. Provides a permanent 24/7 URL (e.g. \https://clairescope.streamlit.app\).
