#!/usr/bin/env python3
"""Mail the submitter of a proteocast.ijm.fr job when it ends (optional email field, user 2026-11-19).

Called as the last line of each job's Slurm script (as www-data): notify_user.py <job id>. The address is in
NOTIFY/<job id> (written by views.handle_upload, outside /data/jobs because Apache serves the job folders; job_alert.py
deletes it after 30 days). No address: nothing to do. Success = a ProteoCast output file in the job folder.
Sender: /data/msa_store/mail_from if present (one line, e.g. "ProteoCast <x@y>"), else Postfix's default.
Standard library only.
"""
import os, subprocess, sys
from email.message import EmailMessage

NOTIFY = '/var/lib/proteocast-notify'
JOBS = '/data/jobs'
FROM_FILE = '/data/msa_store/mail_from'
SITE = 'https://proteocast.ijm.fr/results/results?q=job'


def main():
    j = sys.argv[1] if len(sys.argv) > 1 else ''
    if not (j.isdigit() and len(j) == 14):
        return 1
    p = os.path.join(NOTIFY, j)
    if not os.path.exists(p):
        return 0
    to = open(p).read().strip()
    d = os.path.join(JOBS, j)
    files = os.listdir(d)
    m = EmailMessage()
    m['To'] = to
    if os.path.exists(FROM_FILE):
        m['From'] = open(FROM_FILE).read().strip()
    link = SITE + j
    if any('ProteoCast' in f for f in files):
        m['Subject'] = f'ProteoCast: your job {j} is finished'
        body = f'Your ProteoCast job {j} is finished. The results are here:\n\n  {link}\n\nYou can download all output files from that page.\n'
    else:
        reason = 'The job ended without results.'
        if 'msa_error.txt' in files:
            reason = open(os.path.join(d, 'msa_error.txt'), errors='replace').read().strip()
        m['Subject'] = f'ProteoCast: your job {j} failed'
        body = (f'Your ProteoCast job {j} did not produce results.\n\n  {reason}\n\n  {link}\n\n'
                'The ProteoCast team has been notified. You can submit the job again, or upload your own alignment.\n')
    m.set_content(body + '\n--\nProteoCast web server, https://proteocast.ijm.fr\n'
                  'You receive this message because this address was given when the job was submitted.\n')
    subprocess.run(['/usr/sbin/sendmail', '-t', '-oi'], input=m.as_bytes(), check=True)
    return 0


if __name__ == '__main__':
    sys.exit(main())
