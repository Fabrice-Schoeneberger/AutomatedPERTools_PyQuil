FROM rigetti/forest:latest

WORKDIR /usr/local/bin

COPY . .

CMD ["python", "pyquil_program.py"]