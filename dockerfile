FROM rigetti/forest:latest

WORKDIR /usr/local/bin

COPY . .

CMD ["python", "TrotterExample.py", "--onlyTomography", "--pntsamples", "4", "--pntsinglesamples", "10"]