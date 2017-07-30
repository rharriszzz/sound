/*
gcc -g -O2 -std=gnu99 -Wall -o harmonics.exe harmonics.c -lFLAC -lvorbisenc -lvorbis -logg -L/usr/local/lib -llibsndfile
*/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <ctype.h>

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
  for (int harmonic=1; harmonic<=8; harmonic++) {
    double h = 12 * log((double)harmonic) / log(2);
    int o = (int)((h+0.5)/12);
    printf("%d %d %6.3f\n", harmonic, o, h - o*12);
  }
  for (double note_number = 0; note_number<=88; note_number++) {
    char note_name[13];
    note_number_to_note_name(note_number, note_name, sizeof(note_name), TRUE);
    printf("%6.3f %-10s %12.3f\n", note_number, note_name, note_number_to_frequency(note_number));
  }
  for (double note_number = 0; note_number<=88; note_number++) {
    char note_name[13];
    note_number_to_note_name(note_number, note_name, sizeof(note_name), TRUE);
    printf("%6.3f %-10s ", note_number, note_name);
    for (int harmonic=1; harmonic<=18; harmonic++) {
      double harmonic_note_number = note_number + 12 * log((double)harmonic) / log(2);
      char harmonic_note_name[13];
      note_number_to_note_name(harmonic_note_number, harmonic_note_name, sizeof(harmonic_note_name), TRUE);
      if (harmonic>1 && (harmonic%6)==1) printf("\n                  ");
      printf(" [%d %6.3f %-10s]", harmonic, harmonic_note_number, harmonic_note_name);
    }
    printf("\n");
  }
}

