import sys
from email.message import EmailMessage
from os.path import dirname, realpath
import os
#from app2 import password
import ssl
import smtplib

path = dirname(realpath(__file__))
while "0_data" not in os.listdir(path):
    path = dirname(path)
path = os.path.abspath(path + "/0_data/credentials/lol.txt")
with open(path, 'r') as file:
    email_password = file.read()

email_sender = 'cdng.assistant@gmail.com'


try:
    email_receiver = sys.argv[1] 
except:
    email_receiver = 'elvisavfc65@gmail.com'

try:
    subject = sys.argv[2]
except:
    subject = 'MATLAB Notification'

try:
    body = sys.argv[3]
except:
    body = """
        This email is to notify you that your matlab code has finished executing.
        """

em = EmailMessage()
em['From'] = email_sender
em['To'] = email_receiver
em['Subject'] = subject
em.set_content(body)

context = ssl.create_default_context()

with smtplib.SMTP_SSL('smtp.gmail.com', 465, context = context) as smtp:
    smtp.login(email_sender, email_password)

    smtp.sendmail(email_sender, email_receiver, em.as_string())

