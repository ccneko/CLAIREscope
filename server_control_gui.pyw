import os
import subprocess
import re
import webbrowser
import tkinter as tk
from tkinter import messagebox
from tkinter import ttk
from PIL import Image, ImageDraw, ImageTk

# WSL paths for CLAIREscope
APP_PY_PATH = "/home/claire/dev/CLAIREscope/app.py"
STREAMLIT_BIN = "/home/claire/Software/pyenvs/bioinfo/.venv/bin/streamlit"

def create_rounded_rect_image(width, height, radius, fill_color, stroke_color=None, stroke_width=1):
    image = Image.new("RGBA", (width, height), (0, 0, 0, 0))
    draw = ImageDraw.Draw(image)
    draw.rounded_rectangle(
        (0, 0, width-1, height-1),
        radius=radius,
        fill=fill_color,
        outline=stroke_color,
        width=stroke_width
    )
    return image

class ServerManagerApp:
    def __init__(self, root):
        self.root = root
        self.root.title("CLAIREscope Server Manager")
        self.root.geometry("520x460")
        self.root.resizable(False, False)
        
        # Configure root background (Soft light-gray/slate background)
        self.root.configure(bg="#F8FAFC")
        
        # Dropdown listbox popup customizations
        self.root.option_add('*TCombobox*Listbox.font', ('Segoe UI', 10))
        self.root.option_add('*TCombobox*Listbox.background', '#FFFFFF')
        self.root.option_add('*TCombobox*Listbox.selectBackground', '#0F172A')
        self.root.option_add('*TCombobox*Listbox.selectForeground', '#FFFFFF')
        self.root.option_add('*TCombobox*Listbox.relief', 'flat')
        self.root.option_add('*TCombobox*Listbox.borderWidth', '0')
        
        # Style configuration
        self.style = ttk.Style()
        self.style.theme_use('clam')
        
        # Global background & font config
        self.style.configure('.', font=('Segoe UI', 10), background='#F8FAFC')
        self.style.configure('TFrame', background='#F8FAFC')
        
        # Typography
        self.style.configure('Header.TLabel', font=('Segoe UI', 15, 'bold'), foreground='#0F172A', background='#F8FAFC')
        self.style.configure('Section.TLabel', font=('Segoe UI', 11, 'bold'), foreground='#1E293B', background='#FFFFFF')
        self.style.configure('Normal.TLabel', font=('Segoe UI', 10), foreground='#475569', background='#FFFFFF')
        self.style.configure('Status.TLabel', font=('Segoe UI', 9, 'italic'), foreground='#64748B', background='#F8FAFC')
        
        # Combobox style (Clam theme)
        self.style.configure('TCombobox', fieldbackground='#FFFFFF', bordercolor='#E2E8F0', lightcolor='#E2E8F0', darkcolor='#E2E8F0', arrowcolor='#0F172A', borderwidth=1, relief='solid')

        # Layout Main Frame
        main_frame = ttk.Frame(root, padding="20")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Title Header
        title_label = ttk.Label(main_frame, text="Server Manager", style='Header.TLabel')
        title_label.grid(row=0, column=0, columnspan=3, pady=(0, 15), sticky=tk.W)
        
        # CARD 1: Launch Instance Card (480 x 100, corner radius 8)
        self.card_launch_img = create_rounded_rect_image(480, 100, 8, "#FFFFFF", "#E2E8F0")
        self.card_launch_photo = ImageTk.PhotoImage(self.card_launch_img)
        
        self.card_launch = tk.Label(main_frame, image=self.card_launch_photo, bg="#F8FAFC", borderwidth=0)
        self.card_launch.image = self.card_launch_photo
        self.card_launch.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 15))
        self.card_launch.grid_propagate(False)
        
        # Card 1 Content (placed absolutely)
        lbl_start = ttk.Label(self.card_launch, text="Start Instance", style='Section.TLabel')
        lbl_start.place(x=20, y=12)
        
        lbl_port = ttk.Label(self.card_launch, text="Port:", style='Normal.TLabel')
        lbl_port.place(x=20, y=52)
        
        # Rounded Entry Box (Bootstrap .form-control style: white, CCC border, radius 4)
        self.entry_bg_img = create_rounded_rect_image(120, 32, 4, "#FFFFFF", "#CCCCCC")
        self.entry_bg_photo = ImageTk.PhotoImage(self.entry_bg_img)
        self.entry_container = tk.Label(self.card_launch, image=self.entry_bg_photo, bg="#FFFFFF", borderwidth=0)
        self.entry_container.image = self.entry_bg_photo
        self.entry_container.place(x=70, y=47, width=120, height=32)
        
        self.port_entry = tk.Entry(self.entry_container, bg="#FFFFFF", fg="#555555", font=("Segoe UI", 10), borderwidth=0, highlightthickness=0, insertbackground="#555555")
        self.port_entry.insert(0, "8501")
        self.port_entry.pack(padx=10, pady=6, fill=tk.BOTH, expand=True)
        
        # Rounded Launch Button (White font on Slate)
        self.btn_start = self.make_rounded_button(self.card_launch, "Launch", "#0F172A", "#1E293B", self.start_server, 100, 32)
        self.btn_start.place(x=360, y=47)
        
        # CARD 2: Active Instances Card (480 x 180, corner radius 8)
        self.card_active_img = create_rounded_rect_image(480, 180, 8, "#FFFFFF", "#E2E8F0")
        self.card_active_photo = ImageTk.PhotoImage(self.card_active_img)
        
        self.card_active = tk.Label(main_frame, image=self.card_active_photo, bg="#F8FAFC", borderwidth=0)
        self.card_active.image = self.card_active_photo
        self.card_active.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=(0, 10))
        self.card_active.grid_propagate(False)
        
        # Card 2 Content
        lbl_run = ttk.Label(self.card_active, text="Running WSL Instances", style='Section.TLabel')
        lbl_run.place(x=20, y=12)
        
        # Dropdown / Combobox (spacious with 20px gap before Refresh)
        self.servers_combo = ttk.Combobox(self.card_active, state="readonly")
        self.servers_combo.place(x=20, y=50, width=320, height=32)
        
        # Rounded Refresh Button (White font on Gray)
        self.btn_refresh = self.make_rounded_button(self.card_active, "Refresh", "#64748B", "#475569", self.refresh_servers, 100, 32)
        self.btn_refresh.place(x=360, y=50)
        
        # Rounded Open Browser Button (White font on Slate)
        self.btn_open = self.make_rounded_button(self.card_active, "Open Browser", "#0F172A", "#1E293B", self.open_browser, 140, 34)
        self.btn_open.place(x=20, y=110)
        
        # Rounded Stop Server Button (White font on Red)
        self.btn_stop = self.make_rounded_button(self.card_active, "Stop Server", "#EF4444", "#DC2626", self.stop_server, 120, 34)
        self.btn_stop.place(x=175, y=110)
        
        # Status Bar
        self.status_var = tk.StringVar()
        self.status_var.set("Ready")
        status_bar = ttk.Label(main_frame, textvariable=self.status_var, style='Status.TLabel')
        status_bar.grid(row=3, column=0, columnspan=3, pady=(5, 0), sticky=tk.W)
        
        # Initial population of running servers
        self.running_servers = []
        self.refresh_servers()

    def make_rounded_button(self, parent, text, bg_color, active_bg_color, command, width, height, radius=4):
        normal_img = create_rounded_rect_image(width, height, radius, bg_color)
        active_img = create_rounded_rect_image(width, height, radius, active_bg_color)
        
        photo_normal = ImageTk.PhotoImage(normal_img)
        photo_active = ImageTk.PhotoImage(active_img)
        
        btn = tk.Button(
            parent,
            image=photo_normal,
            text=text,
            compound='center',
            font=('Segoe UI', 9, 'bold'),
            fg='#FFFFFF',
            borderwidth=0,
            highlightthickness=0,
            activebackground="#FFFFFF",
            bg="#FFFFFF",
            command=command,
            cursor='hand2'
        )
        btn.image = photo_normal
        btn.active_image = photo_active
        
        btn.bind("<Enter>", lambda e: btn.config(image=photo_active))
        btn.bind("<Leave>", lambda e: btn.config(image=photo_normal))
        
        return btn

    def get_running_servers(self):
        try:
            creationflags = 0x08000000 if os.name == 'nt' else 0
            output = subprocess.check_output(
                ["wsl", "ps", "aux"],
                stderr=subprocess.DEVNULL,
                creationflags=creationflags
            ).decode("utf-8", errors="ignore")
        except Exception:
            self.status_var.set("Error: Failed to query WSL.")
            return []
        
        servers = []
        for line in output.splitlines():
            if "streamlit run" in line and ("CLAIREscope" in line or "app.py" in line):
                parts = line.split(None, 10)
                if len(parts) >= 11:
                    pid = parts[1]
                    cmd = parts[10]
                    port_match = re.search(r"--server\.port\s+(\d+)", cmd)
                    port = port_match.group(1) if port_match else "8501"
                    servers.append({
                        "pid": pid,
                        "port": port,
                        "cmd": cmd
                    })
        return sorted(servers, key=lambda x: x["port"])

    def refresh_servers(self):
        self.running_servers = self.get_running_servers()
        if not self.running_servers:
            self.servers_combo['values'] = ("No active instances",)
            self.servers_combo.current(0)
            self.status_var.set("No active servers running.")
            self.btn_open.config(state='disabled')
            self.btn_stop.config(state='disabled')
        else:
            options = [f"Port {s['port']} (PID: {s['pid']})" for s in self.running_servers]
            self.servers_combo['values'] = options
            self.servers_combo.current(0)
            self.status_var.set(f"Found {len(self.running_servers)} active server(s).")
            self.btn_open.config(state='normal')
            self.btn_stop.config(state='normal')

    def start_server(self):
        port = self.port_entry.get().strip()
        if not port.isdigit():
            messagebox.showerror("Error", "Port must be a valid integer.")
            return
            
        cmd = [
            "wsl", STREAMLIT_BIN, "run", APP_PY_PATH,
            "--server.port", port,
            "--server.address", "0.0.0.0",
            "--server.headless", "true"
        ]
        
        try:
            creationflags = 0x08000000 if os.name == 'nt' else 0
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                creationflags=creationflags
            )
            if not hasattr(self, 'launched_processes'):
                self.launched_processes = []
            self.launched_processes.append(proc)
            self.status_var.set(f"Launching instance on port {port}...")
            self.root.after(2500, self.refresh_servers)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch: {str(e)}")

    def open_browser(self):
        if self.btn_open.cget('state') == 'disabled':
            return
        selected_idx = self.servers_combo.current()
        if not self.running_servers or selected_idx < 0:
            return
        selected_server = self.running_servers[selected_idx]
        port = selected_server['port']
        url = f"http://localhost:{port}"
        self.status_var.set(f"Opening browser for port {port}...")
        webbrowser.open(url)

    def stop_server(self):
        if self.btn_stop.cget('state') == 'disabled':
            return
        selected_idx = self.servers_combo.current()
        if not self.running_servers or selected_idx < 0:
            return
            
        selected_server = self.running_servers[selected_idx]
        pid = selected_server['pid']
        port = selected_server['port']
        
        try:
            creationflags = 0x08000000 if os.name == 'nt' else 0
            subprocess.call(
                ["wsl", "kill", "-9", pid],
                stderr=subprocess.DEVNULL,
                creationflags=creationflags
            )
            self.status_var.set(f"Terminated instance on port {port} (PID: {pid}).")
            self.root.after(1000, self.refresh_servers)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to stop: {str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = ServerManagerApp(root)
    root.mainloop()
