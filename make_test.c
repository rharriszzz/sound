/*
gcc -g -O2 -std=gnu99 -Wall -o make_test.exe make_test.c -lFLAC -lvorbisenc -lvorbis -logg -L/usr/local/lib -llibsndfile

./make_test.exe -f /cygdrive/c/users/rharris/Music/test1.wav -fr 100000
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <complex.h>
#include <time.h>
#include <sndfile.h>

#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

int note_name_to_note_number(char *note_name);
void note_number_to_note_name(double note_number_double, char *note_name, int note_name_size, int show_remainder_in_cents);
double note_number_to_frequency(double note_number);
double frequency_to_note_number(double frequency);

int note_name_to_note_number(char *note_name)
{
  int note_letter = toupper((int)note_name[0]);
  int is_sharp = '#' == note_name[1];
  unsigned char octave_char = is_sharp ? note_name[2] : note_name[1];
  int note_number = is_sharp + \
    (octave_char - '4') * 12 +
    (note_letter - 'A') * 2 +
    ((note_letter >= 'F') ? -16 : (note_letter >= 'C') ? -13 : 0);
  /* C4 = -9, C#4 = -8, D4 = -7, D#4 = -6, E4 = -5, F4 = -4, F#4 = -3, G4 = -2, G#4 = -1, A4 = 0, A#4 = 1, B4 = 2, C5 = 3, C#5 = 4 */
  return note_number + 49;
}

void note_number_to_note_name(double note_number_double, char *note_name, int note_name_size, int show_remainder_in_cents)
{
  int note_number = lrint(note_number_double);
  int octave = (note_number+8)/12;
  int note = (note_number+8)%12;
  const char *note_names[] = {"C_", "C#", "D_", "D#", "E_", "F_", "F#", "G_", "G#", "A_", "A#", "B_"};
  snprintf(note_name, note_name_size, "%s%d", note_names[note], octave);
  if (show_remainder_in_cents) {
    int el = strlen(note_name);
    double remainder = note_number_double - (double)note_number;
    snprintf(note_name+el, note_name_size-el, "%+07.3f", remainder*100);
  }
}

double note_number_to_frequency(double note_number)
{
  return 440.0 * exp(log(2.0) * (note_number - 49) / 12);
}

double frequency_to_note_number(double frequency)
{
  return log(frequency / 440.0) * 12.0 / log(2.0) + 49;
}


int main(int argc, char **argv)
{
  char *filename = NULL;
  int frames = 1024;
  char *note_name = "A4";
  double note_number = 49;
  int samples_per_second = 44100;
               
  for (int p = 1; p < argc; p++) {
    if (0 == strcmp(argv[p], "-f") || 0 == strcmp(argv[p], "--file")) {
      if (++p < argc) filename = argv[p];
    } else if (0 == strcmp(argv[p], "-fr") || 0 == strcmp(argv[p], "--frames")) {
      if (++p < argc) frames = strtol(argv[p], NULL, 10);
   } else if (0 == strcmp(argv[p], "-n") || 0 == strcmp(argv[p], "--note")) {
      if (++p < argc) {
        note_name = argv[p];
        note_number = (double)note_name_to_note_number(note_name);
      }
    } else if (0 == strcmp(argv[p], "-nn") || 0 == strcmp(argv[p], "--notenumber")) {
      if (++p < argc) {
        note_number = strtod(argv[p], NULL);
        int note_name_size = 13;
        char *note_name = (char *)malloc(note_name_size);
        note_number_to_note_name(note_number, note_name, note_name_size, TRUE);
      }
     } else {
      printf("Option %s was seen, but it is not supported\n", argv[p]);
    }
  }
  if (filename == NULL) {
    printf("Please specify a filename, using -f\n");
    exit(0);
  }

  double note = 440 * exp(log(2.0) * (note_number - 49) / 12);
  printf("note_name=%s, note_number %8.3f, freq=%9.3f\n", note_name, note_number, note);

  double *data = (double *)calloc(frames, sizeof(double));
  for (int f=0; f<frames; f++) {
    double sum = 0.0;
    for (int j=1; j<=8; j++) {
      sum += sin(2*M_PI*note*j*f/samples_per_second) * ((j==1) ? 1.0 : (j==2) ? 0.1 : 0.05);
    }
    data[f] = sum * exp(-f*(4.0/samples_per_second)) / 8;
  }
  
  SF_INFO info_s, *info = &info_s;
  memset(info, 0, sizeof(*info));
  info->samplerate = samples_per_second;
  info->channels = 1;
  if (strstr(filename, ".flac"))
    info->format = SF_FORMAT_FLAC + SF_FORMAT_PCM_16; /* or _24 or _32 */
  if (strstr(filename, ".wav"))
    info->format = SF_FORMAT_WAV + SF_FORMAT_PCM_16; /* or _24 or _32 */
  SNDFILE *out = sf_open(filename, SFM_WRITE, info);
  sf_write_double(out, data, frames);
  sf_close(out);
}
