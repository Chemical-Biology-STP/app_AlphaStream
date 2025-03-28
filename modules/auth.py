import secrets
import sqlite3
import subprocess

import bcrypt
import streamlit as st

DB_FILE = "app_AlphaStream/users.db"


# ✅ Initialize Database with Reset Code Column
def init_db():
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute(
        """
        CREATE TABLE IF NOT EXISTS users (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            email TEXT UNIQUE NOT NULL,
            password TEXT NOT NULL,
            verified INTEGER DEFAULT 0,
            verification_code TEXT,
            reset_code TEXT
        )
    """
    )
    conn.commit()
    conn.close()


# ✅ Hash Passwords Securely
def hash_password(password):
    return bcrypt.hashpw(password.encode(), bcrypt.gensalt()).decode()


# ✅ Verify Passwords
def verify_password(password, hashed_password):
    return bcrypt.checkpw(password.encode(), hashed_password.encode())


# ✅ Generate Secure Tokens
def generate_token():
    return secrets.token_urlsafe(16)


# ✅ Send Email Using `sendmail`
def send_email_sendmail(to_email, subject, body):
    try:
        email_content = f"To: {to_email}\nSubject: {subject}\n\n{body}"
        process = subprocess.Popen(["/usr/sbin/sendmail", "-t"], stdin=subprocess.PIPE)
        process.communicate(email_content.encode())
        return True, None  # ✅ Now it always returns two values
    except Exception as e:
        return False, str(e)


# ✅ Send Email Using `mailx`
def send_email_mailx(to_email, subject, body):
    try:
        command = f'echo "{body}" | mail -s "{subject}" {to_email}'
        result = subprocess.run(command, shell=True, capture_output=True, text=True)
        return result.returncode == 0, (
            result.stderr.strip() if result.returncode != 0 else None
        )
    except Exception as e:
        return False, str(e)


# ✅ Request Password Reset
def request_password_reset(email):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()

    # Check if user exists
    c.execute("SELECT * FROM users WHERE email = ?", (email,))
    user = c.fetchone()
    if not user:
        conn.close()
        return False, "No account found with this email."

    # Generate reset token & save to database
    reset_code = generate_token()
    c.execute("UPDATE users SET reset_code = ? WHERE email = ?", (reset_code, email))
    conn.commit()
    conn.close()

    # ✅ Construct Reset Email
    reset_link = f"http://10.0.208.175:8501/AlphaStream/?reset={reset_code}"
    email_body = f"Click the link below to reset your password:\n\n{reset_link}"

    # ✅ Use sendmail or mailx (Choose one)
    email_sent, error = send_email_sendmail(
        email, "Reset Your AlphaStream Password", email_body
    )
    # email_sent, error = send_email_mailx(email, "Reset Your AlphaStream Password", email_body)

    if email_sent:
        return True, "A password reset email has been sent to your inbox."
    else:
        return False, f"Error sending reset email: {error}"


# ✅ Reset Password
def reset_password(reset_code, new_password):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()

    # Check if the reset code exists
    c.execute("SELECT email FROM users WHERE reset_code = ?", (reset_code,))
    user = c.fetchone()

    if user:
        email = user[0]

        # ✅ Remove reset code before changing password (prevents reusing the link)
        c.execute("UPDATE users SET reset_code = NULL WHERE email = ?", (email,))
        conn.commit()

        # ✅ Hash and update new password
        hashed_password = hash_password(new_password)
        c.execute(
            "UPDATE users SET password = ? WHERE email = ?", (hashed_password, email)
        )
        conn.commit()
        conn.close()
        return True, "Your password has been reset successfully!"

    conn.close()
    return False, "Invalid or expired password reset link."


def login_user(email, password):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()
    c.execute("SELECT password, verified FROM users WHERE email = ?", (email,))
    user = c.fetchone()
    conn.close()

    if user:
        if user[1] == 0:
            return False, "Your account is not verified. Check your email."
        if verify_password(password, user[0]):
            return True, "Login successful!"
    return False, "Invalid email or password."


def register_user(email, password):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()

    # Check if email is already registered
    c.execute("SELECT * FROM users WHERE email = ?", (email,))
    if c.fetchone():
        conn.close()
        return False, "Email already registered. Please log in."

    # Insert new user
    hashed_password = hash_password(password)
    verification_code = generate_token()
    c.execute(
        "INSERT INTO users (email, password, verified, verification_code) VALUES (?, ?, ?, ?)",
        (email, hashed_password, 0, verification_code),
    )
    conn.commit()
    conn.close()

    # Send Verification Email
    email_sent = send_email_sendmail(
        email,
        "Verify Your AlphaStream Account",
        f"Click the link to verify your account:\n\nhttp://10.0.208.175:8501/AlphaStream/?verify={verification_code}",
    )
    if email_sent:
        return True, "Registration successful! Check your email to verify your account."
    else:
        return False, "Error sending verification email."


def verify_user(verification_code):
    conn = sqlite3.connect(DB_FILE)
    c = conn.cursor()

    # Check if the verification code exists
    c.execute(
        "SELECT email, verified FROM users WHERE verification_code = ?",
        (verification_code,),
    )
    user = c.fetchone()

    if user:
        email, verified = user
        if verified == 1:
            conn.close()
            return False, "This verification link has already been used."

        # ✅ Remove verification code and mark the user as verified
        c.execute(
            "UPDATE users SET verified = 1, verification_code = NULL WHERE email = ?",
            (email,),
        )
        conn.commit()
        conn.close()
        return True, "Your account has been successfully verified!"

    conn.close()
    return False, "Invalid or expired verification link."


# ✅ Initialize the database
init_db()
