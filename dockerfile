FROM rigetti/forest:latest

WORKDIR /usr/local/bin

COPY . .

CMD ["python", "test.py"]