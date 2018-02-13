import files;

app (file o) hostname() {
  "hostname" @stdout=o;
}

trace(filename(hostname()));
trace(read(hostname()));
