import sqlite3  # ✅ Ensure sqlite3 is imported
import time

import streamlit as st

from modules.auth import (
    login_user,
    register_user,
    request_password_reset,
    reset_password,
    verify_user,
)
from modules.ui_components.tabs import render_tabs


# ✅ Use JavaScript for same-tab redirection
def redirect_same_tab(url):
    """Redirect within the same tab using JavaScript."""
    st.markdown(
        f"""
        <meta http-equiv="refresh" content="0; URL={url}">
        """,
        unsafe_allow_html=True,
    )


def main():
    """
    Streamlit app entry point. Handles user authentication (registration, login, logout, password reset).
    """
    st.title("AlphaStream")

    # ✅ Show Latest Changelog (Expanded by Default)
    with st.expander("📜 Latest Changelog (7 February 2024)", expanded=False):
        st.markdown(
            """
            ### 🔹 Latest Updates – 7 February 2024  

            #### **1️⃣ View Amino Acid Positions in a Structured Format**  
            - You can now **see the position of each amino acid in your input sequence** in a structured table.  
            - The table is **organized in rows of 10 amino acids**, making it easier to reference specific positions.  
            - This feature is available under a **collapsible section** labeled **“📍 Show Sequence with Positions”**.  
            - Click on the section to expand it and view the table.  

            **📌 Why this helps:**  
            🔹 Makes it easier to **select the correct position** for PTM modifications.  
            🔹 Ensures that **each amino acid position is clear and structured**.  

            ---

            #### **2️⃣ Easily Select PTMs Based on Amino Acids**  
            - Instead of manually entering a PTM **CCD modification code**, you can now:  
            ✅ **Select an amino acid from a dropdown menu**.  
            ✅ **Choose from a list of valid PTM modifications** specific to that amino acid.  
            - The **dropdown menu is now sorted alphabetically**, making it **easier to find an amino acid**.  

            **📌 Why this helps:**  
            🔹 No need to **remember or look up modification codes** manually.  
            🔹 Faster and more accurate **PTM selection**.  
            🔹 **No more scrolling through a long, unordered list**—everything is neatly sorted!  

            ---

            #### **3️⃣ Clickable Links for PTM Information**  
            - Once a PTM modification is selected, a **clickable link** appears next to it.  
            - This **redirects you to a trusted database** where you can **view details about the PTM**.  
            - Example: Selecting **HY3** will generate a link like:  
            👉 **[🔗 More Info on HY3](https://www.ebi.ac.uk/pdbe-srv/pdbechem/chemicalCompound/show/HY3)**  

            **📌 Why this helps:**  
            🔹 Provides **official information** about each PTM.  
            🔹 No need to search manually—just click and learn!  

            ---

            ### **🔹 Summary of Benefits**  
            ✅ **Amino acid positions are now clearly displayed** for better accuracy.  
            ✅ **PTM selection is faster and more user-friendly**, with an **alphabetically sorted amino acid list**.  
            ✅ **Clickable PTM links make it easy to access more information** about each modification.  

            These improvements were made based on **user feedback**, and we hope they make your experience **simpler, faster, and more efficient**! 🚀  

            If you have any further suggestions, feel free to share them with us! 😊
            """
        )

    # ✅ Show Past Changelogs (Collapsed by Default)
    with st.expander("📜 Past Changelogs"):
        st.markdown(
            """
            ### 🔹 Changelog – 6 February 2024  

            #### **🔐 Implementation of User Accounts for Enhanced Security & Privacy**  
            - Users must now **create an account with a password**, ensuring only the rightful owner of an email can log in.  
            - Email verification is **required** before an account can be used, **preventing unauthorized access**.  
            - A **password reset feature** has been added, allowing users to securely regain access to their accounts.  

            ---

            #### **🚀 Improved Job Submission & Status Checking**  
            *Credits to: Celine Bouchoux ([Celine.Bouchoux@crick.ac.uk](mailto:Celine.Bouchoux@crick.ac.uk)) and Dave Briggs ([david.briggs@crick.ac.uk](mailto:david.briggs@crick.ac.uk))*  
            - Submitting jobs, checking job statuses, and retrieving files is now **faster and more efficient**.  

            ---

            #### **✅ Ensuring Valid Job Names**  
            *Credits to: Jack Dainton ([jack.dainton@crick.ac.uk](mailto:jack.dainton@crick.ac.uk))*  
            - Job names now undergo validation to **prevent errors and ensure correct formatting**.  

            ---

            #### **🔄 Users Can Now Check Individual Job Statuses or Refresh All Job Statuses at Once**  
            - This allows for **more flexibility when tracking job progress**.  

            ---

            #### **🔐 Enhanced Email Verification Redirect for Improved Security and Privacy**  
            - Users are now redirected to `http://10.0.208.175:8501/AlphaStream` **without the verification link remaining in the address bar**.  
            - This **prevents unauthorized access to verification links** and **protects user privacy**.  
            - Users no longer see the **"Invalid verification link or already verified"** message when refreshing the page.  

            ---

            #### **🔒 Password Reset Links Now Expire Immediately After Use**  
            - Previously, users could **reuse old password reset links**, which was a security risk.  
            - Now, once a password is reset, **the reset link becomes invalid immediately**.  
            - Users attempting to reuse an expired link will see **"Invalid or expired password reset link."**  

            ---

            #### **🛠️ Redirects Now Occur in the Same Tab**  
            - Previously, verification and password reset links opened in a **new tab**.  
            - Now, all redirects occur within the **same tab for a smoother user experience**.  
            """
        )

    # ✅ Ensure session state is initialized
    if "logged_in" not in st.session_state:
        st.session_state.logged_in = False
        st.session_state.email = ""

    # ✅ Handle Password Reset Links (Show Reset Page Only)
    reset_code = st.query_params.get("reset")
    if reset_code:
        st.subheader("🔑 Reset Your Password")

        # ✅ Check if reset code is still valid before showing the form
        conn = sqlite3.connect("users.db")
        c = conn.cursor()
        c.execute("SELECT email FROM users WHERE reset_code = ?", (reset_code,))
        valid_code = c.fetchone()
        conn.close()

        if not valid_code:
            st.error("❌ Invalid or expired password reset link.")
            return  # ✅ Stop execution here

        # ✅ Display reset password form
        new_password = st.text_input("Enter new password:", type="password")
        confirm_password = st.text_input("Confirm new password:", type="password")

        if st.button("Reset Password"):
            if new_password != confirm_password:
                st.error("❌ Passwords do not match.")
            else:
                success, message = reset_password(reset_code, new_password)
                if success:
                    st.success(
                        "✅ Password reset successful! Redirecting to login page..."
                    )

                    # ✅ Wait for 3 seconds before redirecting
                    time.sleep(3)

                    # ✅ Redirect to login page in the same tab
                    st.query_params.clear()  # Clear the URL query string
                    redirect_same_tab(
                        "http://10.0.208.175:8501/AlphaStream"
                    )  # ✅ Redirect in the same tab
                    st.rerun()
                else:
                    st.error(message)

        return  # ✅ Prevent login/register options from appearing

    # ✅ Handle Email Verification Links
    verification_code = st.query_params.get("verify")
    if verification_code:
        success, message = verify_user(verification_code)
        if success:
            st.success(message)
            st.info("🔄 Redirecting to the login page...")

            # ✅ Wait 3 seconds before redirecting
            time.sleep(3)

            # ✅ Redirect to login page in the same tab
            st.query_params.clear()  # Clear the query string
            redirect_same_tab(
                "http://10.0.208.175:8501/AlphaStream"
            )  # ✅ Redirect in the same tab
            st.rerun()
        else:
            st.error(message)

    # ✅ If the user is logged in, render the tabs
    if st.session_state.logged_in:
        st.success(f"✅ Logged in as {st.session_state.email}")
        render_tabs(st.session_state.email)

        # ✅ Logout Button
        if st.button("Logout"):
            st.session_state.logged_in = False
            st.session_state.email = ""
            st.rerun()

        return  # ✅ Prevent login form from appearing when logged in

    # ✅ Show registration, login, or forgot password options
    option = st.selectbox("Select an option:", ["Login", "Register", "Forgot Password"])

    if option == "Register":
        st.subheader("🔐 Register")
        email = st.text_input("Enter your email:")
        password = st.text_input("Enter a password:", type="password")
        confirm_password = st.text_input("Confirm password:", type="password")

        if st.button("Register"):
            if not email.endswith("@crick.ac.uk"):
                st.error(
                    "❌ Registration is only allowed for Crick email addresses (`@crick.ac.uk`)."
                )
            elif password != confirm_password:
                st.error("❌ Passwords do not match.")
            elif email and password:
                success, message = register_user(email, password)
                if success:
                    st.success(message)
                    st.info(
                        "📧 Please check your **Inbox, Junk, or Spam** folder for the verification email."
                    )
                else:
                    st.error(message)

    elif option == "Forgot Password":
        st.subheader("🔑 Forgot Password")
        email = st.text_input("Enter your email to reset password:")

        if st.button("Send Reset Link"):
            if not email.endswith("@crick.ac.uk"):
                st.error(
                    "❌ Password reset is only available for Crick email addresses (`@crick.ac.uk`)."
                )
            else:
                success, message = request_password_reset(email)
                if success:
                    st.success(message)
                    st.info(
                        "📧 Please check your **Inbox, Junk, or Spam** folder for the password reset email."
                    )
                else:
                    st.error(message)

    elif option == "Login":
        st.subheader("🔑 Login")

        # ✅ Ensure input fields are rendered
        email = st.text_input("Enter your email:").lower()
        password = st.text_input("Enter your password:", type="password")

        if st.button("Login"):
            if not email.endswith("@crick.ac.uk"):
                st.error(
                    "❌ Login is only allowed for Crick email addresses (`@crick.ac.uk`)."
                )
            else:
                success, message = login_user(email, password)
                if success:
                    st.success(message)
                    st.session_state.logged_in = True
                    st.session_state.email = email
                    st.rerun()  # ✅ Ensure state updates properly
                else:
                    st.error(message)


if __name__ == "__main__":
    main()
