"""
CLAIREscope Desktop Server Controller GUI
Provides a clean, cross-platform Tkinter interface to start, stop, and manage the Streamlit server.
"""

import os
import sys
import time
import signal
import socket
import webbrowser
import subprocess
import threading
import collections
from typing import Optional

import tkinter as tk
from tkinter import ttk, messagebox
from tkinter.scrolledtext import ScrolledText

APP_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP_PY_PATH = os.path.join(APP_DIR, "app.py")


def is_port_in_use(port: int, host: str = "127.0.0.1") -> bool:
    """Check if a TCP port is currently open and bound."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.5)
        return s.connect_ex((host, port)) == 0


class ServerControllerGUI:
    def __init__(self, root: tk.Tk):
        self.root = root
        self.root.title("CLAIREscope — Server Controller")
        self.root.minsize(480, 260)
        self.root.resizable(True, True)

        self.process: Optional[subprocess.Popen] = None
        self.log_history = collections.deque(maxlen=2000)
        self.is_starting = False
        self.log_reader_thread: Optional[threading.Thread] = None

        # Variables
        self.port_var = tk.StringVar(value="8501")
        self.listen_lan_var = tk.BooleanVar(value=True)
        self.auto_open_var = tk.BooleanVar(value=True)
        self.debug_mode_var = tk.BooleanVar(value=False)
        self.status_text_var = tk.StringVar(value="🔴 Server Stopped")
        self.url_text_var = tk.StringVar(value="http://localhost:8501")

        self._setup_style()
        self._build_ui()

        # Handle window close
        self.root.protocol("WM_DELETE_WINDOW", self.on_close)

        # Periodic status checker
        self.root.after(800, self._periodic_check)

    def _setup_style(self):
        style = ttk.Style()
        try:
            style.theme_use("clam")
        except Exception:
            pass

        style.configure("Title.TLabel", font=("Segoe UI", 14, "bold"), foreground="#1F2937")
        style.configure("SubTitle.TLabel", font=("Segoe UI", 9), foreground="#6B7280")
        style.configure("Status.TLabel", font=("Segoe UI", 11, "bold"))
        style.configure("URL.TLabel", font=("Segoe UI", 10, "underline"), foreground="#2563EB")
        style.configure("Accent.TButton", font=("Segoe UI", 10, "bold"))

    def _build_ui(self):
        # Main container with padding
        self.main_frame = ttk.Frame(self.root, padding="16 16 16 16")
        self.main_frame.pack(fill=tk.BOTH, expand=True)

        # Header
        header_frame = ttk.Frame(self.main_frame)
        header_frame.pack(fill=tk.X, pady=(0, 10))

        title_lbl = ttk.Label(header_frame, text="🔬 CLAIREscope Server Controller", style="Title.TLabel")
        title_lbl.pack(anchor=tk.W)

        subtitle_lbl = ttk.Label(
            header_frame,
            text="Single-Cell & Spatial Transcriptomics Analysis Platform",
            style="SubTitle.TLabel"
        )
        subtitle_lbl.pack(anchor=tk.W)

        # Status Card Box
        self.status_card = tk.LabelFrame(
            self.main_frame,
            text=" Server Status ",
            font=("Segoe UI", 9, "bold"),
            padx=12,
            pady=10,
            bg="#F9FAFB",
            fg="#374151"
        )
        self.status_card.pack(fill=tk.X, pady=(0, 12))

        c_status_row = tk.Frame(self.status_card, bg="#F9FAFB")
        c_status_row.pack(fill=tk.X)

        self.status_badge = tk.Label(
            c_status_row,
            textvariable=self.status_text_var,
            font=("Segoe UI", 11, "bold"),
            bg="#FEE2E2",
            fg="#991B1B",
            padx=10,
            pady=4,
            relief=tk.FLAT
        )
        self.status_badge.pack(side=tk.LEFT)

        self.url_label = tk.Label(
            c_status_row,
            textvariable=self.url_text_var,
            font=("Segoe UI", 10, "underline"),
            fg="#2563EB",
            bg="#F9FAFB",
            cursor="hand2",
            padx=10
        )
        self.url_label.pack(side=tk.LEFT, padx=10)
        self.url_label.bind("<Button-1>", lambda e: self.open_browser())

        # Action Buttons
        btn_frame = ttk.Frame(self.main_frame)
        btn_frame.pack(fill=tk.X, pady=(0, 12))

        self.btn_start = tk.Button(
            btn_frame,
            text="▶ Start Server",
            font=("Segoe UI", 10, "bold"),
            bg="#10B981",
            fg="white",
            activebackground="#059669",
            activeforeground="white",
            relief=tk.FLAT,
            padx=14,
            pady=6,
            command=self.start_server
        )
        self.btn_start.pack(side=tk.LEFT, padx=(0, 6))

        self.btn_stop = tk.Button(
            btn_frame,
            text="⏹ Stop Server",
            font=("Segoe UI", 10, "bold"),
            bg="#EF4444",
            fg="white",
            activebackground="#DC2626",
            activeforeground="white",
            relief=tk.FLAT,
            padx=14,
            pady=6,
            state=tk.DISABLED,
            command=self.stop_server
        )
        self.btn_stop.pack(side=tk.LEFT, padx=(0, 6))

        self.btn_restart = tk.Button(
            btn_frame,
            text="🔄 Restart",
            font=("Segoe UI", 9),
            bg="#F3F4F6",
            fg="#374151",
            activebackground="#E5E7EB",
            relief=tk.FLAT,
            padx=10,
            pady=6,
            state=tk.DISABLED,
            command=self.restart_server
        )
        self.btn_restart.pack(side=tk.LEFT, padx=(0, 6))

        self.btn_browser = tk.Button(
            btn_frame,
            text="🌐 Open in Browser",
            font=("Segoe UI", 9),
            bg="#3B82F6",
            fg="white",
            activebackground="#2563EB",
            activeforeground="white",
            relief=tk.FLAT,
            padx=12,
            pady=6,
            state=tk.DISABLED,
            command=self.open_browser
        )
        self.btn_browser.pack(side=tk.LEFT, padx=(0, 6))

        # Options / Configuration Row
        opt_frame = ttk.Frame(self.main_frame)
        opt_frame.pack(fill=tk.X, pady=(0, 8))

        ttk.Label(opt_frame, text="Port:").pack(side=tk.LEFT, padx=(0, 4))
        self.port_entry = ttk.Entry(opt_frame, textvariable=self.port_var, width=6)
        self.port_entry.pack(side=tk.LEFT, padx=(0, 12))

        self.chk_lan = ttk.Checkbutton(
            opt_frame,
            text="Allow LAN / Remote (0.0.0.0)",
            variable=self.listen_lan_var
        )
        self.chk_lan.pack(side=tk.LEFT, padx=(0, 12))

        self.chk_auto_open = ttk.Checkbutton(
            opt_frame,
            text="Auto-open browser",
            variable=self.auto_open_var
        )
        self.chk_auto_open.pack(side=tk.LEFT, padx=(0, 12))

        # Debug Mode Toggle (Defaults to False / Unchecked)
        self.chk_debug = ttk.Checkbutton(
            opt_frame,
            text="🐞 Debug Mode",
            variable=self.debug_mode_var,
            command=self._toggle_debug_console
        )
        self.chk_debug.pack(side=tk.RIGHT)

        # Collapsible Debug Log Console Container (Hidden by default)
        self.log_container = ttk.Frame(self.main_frame)
        # Not packed initially since debug_mode_var is False

        log_header = ttk.Frame(self.log_container)
        log_header.pack(fill=tk.X, pady=(4, 2))
        ttk.Label(log_header, text="Real-time Server Logs:", font=("Segoe UI", 8, "bold"), foreground="#4B5563").pack(side=tk.LEFT)
        btn_clear = ttk.Button(log_header, text="Clear", width=6, command=self.clear_logs)
        btn_clear.pack(side=tk.RIGHT)

        self.log_text = ScrolledText(
            self.log_container,
            height=10,
            bg="#111827",
            fg="#10B981",
            insertbackground="white",
            font=("Consolas", 8),
            wrap=tk.WORD
        )
        self.log_text.pack(fill=tk.BOTH, expand=True)

    def _toggle_debug_console(self):
        """Show or hide the real-time log console based on Debug Mode checkbox."""
        if self.debug_mode_var.get():
            self.log_container.pack(fill=tk.BOTH, expand=True, pady=(6, 0))
            self._render_all_logs()
        else:
            self.log_container.pack_forget()

    def _render_all_logs(self):
        """Refresh log console with accumulated history."""
        self.log_text.delete("1.0", tk.END)
        for line in self.log_history:
            self.log_text.insert(tk.END, line)
        self.log_text.see(tk.END)

    def _append_log(self, text: str):
        """Append log line to memory history and update text widget if visible."""
        self.log_history.append(text)
        if self.debug_mode_var.get():
            self.log_text.insert(tk.END, text)
            self.log_text.see(tk.END)

    def clear_logs(self):
        self.log_history.clear()
        self.log_text.delete("1.0", tk.END)

    def get_server_url(self) -> str:
        port = self.port_var.get().strip() or "8501"
        return f"http://localhost:{port}"

    def start_server(self):
        """Launch the Streamlit process."""
        if self.process is not None and self.process.poll() is None:
            messagebox.showinfo("Already Running", "CLAIREscope server is already active!")
            return

        port_str = self.port_var.get().strip() or "8501"
        try:
            port = int(port_str)
        except ValueError:
            messagebox.showerror("Invalid Port", "Port must be a valid number (e.g. 8501).")
            return

        if not os.path.exists(APP_PY_PATH):
            messagebox.showerror("File Not Found", f"Cannot find app entrypoint at: {APP_PY_PATH}")
            return

        # Prepare command
        host_ip = "0.0.0.0" if self.listen_lan_var.get() else "127.0.0.1"
        cmd = [
            sys.executable,
            "-m",
            "streamlit",
            "run",
            APP_PY_PATH,
            f"--server.port={port}",
            f"--server.address={host_ip}",
            "--server.headless=true"
        ]

        self.is_starting = True
        self._set_ui_state_starting()
        self._append_log(f"[{time.strftime('%H:%M:%S')}] Starting CLAIREscope on port {port}...\n")
        self._append_log(f"Command: {' '.join(cmd)}\n")

        try:
            # Launch subprocess
            self.process = subprocess.Popen(
                cmd,
                cwd=APP_DIR,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True,
                creationflags=subprocess.CREATE_NEW_PROCESS_GROUP if os.name == 'nt' else 0
            )

            # Start background reader thread
            self.log_reader_thread = threading.Thread(target=self._stream_logs, daemon=True)
            self.log_reader_thread.start()

            # Schedule ready check in 2 seconds
            self.root.after(2000, self._check_server_ready)

        except Exception as e:
            self.is_starting = False
            self.process = None
            self._set_ui_state_stopped()
            self._append_log(f"[ERROR] Failed to launch server: {e}\n")
            messagebox.showerror("Launch Error", f"Failed to start server:\n{e}")

    def _stream_logs(self):
        """Read stdout/stderr lines continuously."""
        if not self.process or not self.process.stdout:
            return
        try:
            for line in iter(self.process.stdout.readline, ''):
                if not line:
                    break
                self.root.after(0, self._append_log, line)
        except Exception:
            pass

    def _check_server_ready(self):
        """Check if server port has become active."""
        if not self.process or self.process.poll() is not None:
            self.is_starting = False
            self._set_ui_state_stopped()
            return

        port = int(self.port_var.get().strip() or "8501")
        if is_port_in_use(port):
            self.is_starting = False
            self._set_ui_state_running()
            if self.auto_open_var.get():
                self.open_browser()
        else:
            # Poll again in 1 second
            if self.is_starting:
                self.root.after(1000, self._check_server_ready)

    def stop_server(self):
        """Terminate the running Streamlit process."""
        if self.process is None:
            self._set_ui_state_stopped()
            return

        self._append_log(f"[{time.strftime('%H:%M:%S')}] Stopping server (PID: {self.process.pid})...\n")
        try:
            if os.name == 'nt':
                # Force kill on Windows process tree
                subprocess.call(['taskkill', '/F', '/T', '/PID', str(self.process.pid)], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            else:
                self.process.send_signal(signal.SIGTERM)
                self.process.wait(timeout=2)
        except Exception:
            try:
                self.process.kill()
            except Exception:
                pass

        self.process = None
        self.is_starting = False
        self._set_ui_state_stopped()
        self._append_log(f"[{time.strftime('%H:%M:%S')}] Server stopped.\n")

    def restart_server(self):
        """Stop then start server."""
        self.stop_server()
        self.root.after(1000, self.start_server)

    def open_browser(self):
        """Open web browser at the active server URL."""
        url = self.get_server_url()
        webbrowser.open(url)

    def _set_ui_state_stopped(self):
        self.status_text_var.set("🔴 Server Stopped")
        self.status_badge.config(bg="#FEE2E2", fg="#991B1B")
        self.url_text_var.set(self.get_server_url())
        self.btn_start.config(state=tk.NORMAL, bg="#10B981")
        self.btn_stop.config(state=tk.DISABLED, bg="#9CA3AF")
        self.btn_restart.config(state=tk.DISABLED)
        self.btn_browser.config(state=tk.DISABLED, bg="#9CA3AF")
        self.port_entry.config(state=tk.NORMAL)

    def _set_ui_state_starting(self):
        self.status_text_var.set("🟡 Starting Server...")
        self.status_badge.config(bg="#FEF3C7", fg="#92400E")
        self.btn_start.config(state=tk.DISABLED, bg="#9CA3AF")
        self.btn_stop.config(state=tk.NORMAL, bg="#EF4444")
        self.btn_restart.config(state=tk.DISABLED)
        self.btn_browser.config(state=tk.DISABLED, bg="#9CA3AF")
        self.port_entry.config(state=tk.DISABLED)

    def _set_ui_state_running(self):
        self.status_text_var.set("🟢 Server Running")
        self.status_badge.config(bg="#D1FAE5", fg="#065F46")
        self.url_text_var.set(self.get_server_url())
        self.btn_start.config(state=tk.DISABLED, bg="#9CA3AF")
        self.btn_stop.config(state=tk.NORMAL, bg="#EF4444")
        self.btn_restart.config(state=tk.NORMAL)
        self.btn_browser.config(state=tk.NORMAL, bg="#3B82F6")
        self.port_entry.config(state=tk.DISABLED)

    def _periodic_check(self):
        """Check if process crashed or exited externally."""
        if self.process is not None and not self.is_starting:
            if self.process.poll() is not None:
                self.process = None
                self._set_ui_state_stopped()
                self._append_log(f"[{time.strftime('%H:%M:%S')}] Server process terminated.\n")
        self.root.after(1000, self._periodic_check)

    def on_close(self):
        """Prompt and cleanup before closing window."""
        if self.process is not None and self.process.poll() is None:
            if messagebox.askyesno("Quit CLAIREscope", "The CLAIREscope server is still running.\n\nDo you want to stop the server and quit?"):
                self.stop_server()
                self.root.destroy()
        else:
            self.root.destroy()


def main():
    root = tk.Tk()
    app = ServerControllerGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
