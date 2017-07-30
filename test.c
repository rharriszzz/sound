/*
gcc -g -O2 -std=gnu99 -Wall -o test.exe test.c fft.c -lFLAC -lvorbisenc -lvorbis -logg -L/usr/local/lib -llibsndfile

./test.exe -f /cygdrive/c/users/rharris/Music/piano-2015-06-14.flac -fft > data_fft/piano-2015-06-14.out

./test.exe -f /cygdrive/c/users/rharris/Music/piano-2015-06-14.flac -i 48 > piano-2015-06-14__48.out
./test.exe -f /cygdrive/c/users/rharris/Music/piano-2015-06-14.flac -cf4 false -f1b 100 > piano-2015-06-14__f1_all.out
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

#include "fft.h"

int use_transform1 = FALSE;
int use_fft = FALSE;
int use_chirp_z = FALSE;
char *directory = "data";
int use_formatted_show_change = FALSE;

int note_name_to_note_number(char *note_name);
void note_number_to_note_name(double note_number_double, char *note_name, int note_name_size, int show_remainder_in_cents);
double note_number_to_frequency(double note_number);
double frequency_to_note_number(double frequency);

long long transform1_calls;
long long transform1_frames;
clock_t transform1_clock;
void transform1(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		double freq,
		double *xkamp_ptr, double *xkphase_ptr, double *xkdb_ptr);

void transformn(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		double note_base, double note_increment, int note_count, int *result_note_count,
		double **note_number_array_ptr, double **frequency_array_ptr,
		double **xkamp_array_ptr, double **xkphase_array_ptr, double **xkdb_array_ptr);

int total_interval_count;
void format_0_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
                       int samples_per_sum, double begin_cutoff, double end_cutoff,
                       int interval_number,
                       int format_0_verbose,
                       int format_1_bins, int fundemental_db,
                       int call_format_1_analysis, int format_1_verbose,
                       int format_4_bins, double format_4_bin_size,
                       int call_format_4_analysis, int format_4_verbose);

void format_1_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, int fundemental_db, int interval_number,
		       double *note_number_ptr, int verbose);

void format_2_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double);

void format_3_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double);

void format_4_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double, int verbose);


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

long long transform1_calls = 0;
long long transform1_frames = 0;
clock_t transform1_clock = 0;

/*

cos(n*x) == 2*cos(x)*cos((n-1)*x) - cos((n-2)*x)
sin(n*x) == 2*cos(x)*sin((n-1)*x) - sin((n-2)*x)

 */
void transform1(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		double freq,
		double *xkamp_ptr, double *xkphase_ptr, double *xkdb_ptr)
{
  clock_t start_clock = clock();
  if (last_frame == 0) last_frame = info->frames - 1;
  int frame_count = 1 + last_frame - first_frame;
  double xkre = 0.0;
  double xkim = 0.0;
  double delta_angle = (2*M_PI*freq/info->samplerate);
  double angle = first_frame * delta_angle;
  int pos = first_frame*info->channels + channel;
  for (int sample=first_frame; sample<=last_frame; sample++) {
    double value = data[pos];
    pos += info->channels;
    double sin_angle, cos_angle;
#if 0
    sin_angle = sin(angle);
    cos_angle = cos(angle);
#else
    sincos(angle, &sin_angle, &cos_angle);
#endif
    xkre += value*cos_angle;
    xkim += value*sin_angle;
    angle += delta_angle;
  }
  if (xkamp_ptr || xkdb_ptr) {
    double xkamp = hypot(xkre, xkim)/frame_count;
    if (xkamp_ptr) {
      *xkamp_ptr = xkamp;
    }
    if (xkdb_ptr) {
      double xkdb = 20 * log10(xkamp);
      *xkdb_ptr = xkdb;
    }
  }
  if (xkphase_ptr) {
    double xkphase = atan2(xkim, xkre);
    *xkphase_ptr = xkphase;
  }
  clock_t end_clock = clock();
  transform1_calls++;
  transform1_frames += frame_count;
  transform1_clock += (end_clock-start_clock);
}

double *maybe_create_double_array(double **array_ptr, int count)
{
  if (0 == array_ptr) return NULL;
  double *array = (double *)calloc(count, sizeof(double));
  *array_ptr = array;
  return array;
}

void transformn(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		double note_base, double note_increment, int note_count,
                int *result_count,
		double **note_number_array_ptr, double **frequency_array_ptr,
		double **xkamp_array_ptr, double **xkphase_array_ptr, double **xkdb_array_ptr)
{
  if (use_transform1) {
    if (result_count) *result_count = note_count;
    double *note_number_array = maybe_create_double_array(note_number_array_ptr, *result_count);
    double *frequency_array = maybe_create_double_array(frequency_array_ptr, *result_count);
    double *xkamp_array = maybe_create_double_array(xkamp_array_ptr, *result_count);
    double *xkphase_array = maybe_create_double_array(xkphase_array_ptr, *result_count);
    double *xkdb_array = maybe_create_double_array(xkdb_array_ptr, *result_count);
    for (int n = 0; n<note_count; n++) {
      double note_number = note_base + n * note_increment;
      if (note_number_array) note_number_array[n] = note_number;
      double frequency = 440.0 * exp(log(2.0) * (note_number - 49) / 12);
      if (frequency_array) frequency_array[n] = frequency;
      transform1(info, data, channel, first_frame, last_frame,
                 frequency,
                 xkamp_array ? &xkamp_array[n] : NULL,
                 xkphase_array ? &xkphase_array[n] : NULL,
                 xkdb_array ? &xkdb_array[n] : NULL);
    }
  } else if (use_fft) {
    /* 1.0 - 0.5, 1.0/100, 88*100; ni*nc = 88 */
    /* harmonic_note_number - 1.5/2, 1.5/600, 600; ni*nc = 1.5 */
    int N = last_frame - first_frame + 1;
    if (result_count) *result_count = N;
    int M = 1;
    while (M < N) {M *= 2;}
    if (N < M) N = M/2;
    int NR = (N+1)/2;
    if (result_count) *result_count = NR;
    complex double *x = (complex double *)calloc(N, sizeof(complex double));
    double *note_number_array = maybe_create_double_array(note_number_array_ptr, NR);
    double *frequency_array = maybe_create_double_array(frequency_array_ptr, NR);
    double *xkamp_array = maybe_create_double_array(xkamp_array_ptr, NR);
    double *xkphase_array = maybe_create_double_array(xkphase_array_ptr, NR);
    double *xkdb_array = maybe_create_double_array(xkdb_array_ptr, NR);
    int pos = first_frame*info->channels + channel;
    int i = 0;
    for (int sample=first_frame; sample<=last_frame; sample++) {
      if (i >= N) break;
      x[i++] = data[pos];
      pos += info->channels;
    }
    fft(x, N, FALSE, FALSE); /* void fft(complex DOUBLE buf[], int n, int inverse, int balanced); */
    double seconds = ((double)N)/info->samplerate;
    for (int i = 0; i < N/2; i++) {
      double frequency = i/seconds;
      if (frequency_array) frequency_array[i] = frequency;
      if (note_number_array) note_number_array[i] = frequency_to_note_number(frequency);
      double xkamp = cabs(x[i]);
      if (xkamp_array) xkamp_array[i] = xkamp;
      if (xkphase_array) xkphase_array[i] = carg(x[i]) / (2 * M_PI);
      if (xkdb_array) xkdb_array[i] = 20 * log10(xkamp);
    }
    free(x);
  } else if (use_chirp_z) {
    /* 1.0 - 0.5, 1.0/100, 88*100; ni*nc = 88 */
    /* harmonic_note_number - 1.5/2, 1.5/600, 600; ni*nc = 1.5 */
    int N = last_frame - first_frame + 1;
    double low_frequency = note_number_to_frequency(note_base);
    double frequency_increment = note_number_to_frequency(note_base+note_increment) - low_frequency;
    double high_frequency = note_number_to_frequency(note_base+note_count*note_increment);
    if (note_count*note_increment > 12) low_frequency = low_frequency / 2;

    int M = high_frequency / frequency_increment;
    double *note_number_array = maybe_create_double_array(note_number_array_ptr, M);
    double *frequency_array = maybe_create_double_array(frequency_array_ptr, M);
    double *xkamp_array = maybe_create_double_array(xkamp_array_ptr, M);
    double *xkphase_array = maybe_create_double_array(xkphase_array_ptr, M);
    double *xkdb_array = maybe_create_double_array(xkdb_array_ptr, M);

    complex double *x = (complex double *)calloc(N, sizeof(complex double));
    complex double *X = (complex double *)calloc(M, sizeof(complex double));
    int pos = first_frame*info->channels + channel;
    int i = 0;
    for (int sample=first_frame; sample<=last_frame; sample++) {
      x[i++] = data[pos];
      pos += info->channels;
    }
    if (result_count) *result_count = N;
    complex double A = 1.0 * cexp(- I * 2 * M_PI * low_frequency / info->samplerate);
    complex double W = 1.0 * cexp(- I * 2 * M_PI * frequency_increment / info->samplerate);
    chirp_z(x, N, A, W, X, M);
    for (int i = 0; i < M; i++) {
      double frequency = low_frequency+i*frequency_increment;
      if (frequency_array) frequency_array[i] = frequency;
      if (note_number_array) note_number_array[i] = frequency_to_note_number(frequency);
      double xkamp = cabs(X[i]);
      if (xkamp_array) xkamp_array[i] = xkamp;
      if (xkphase_array) xkphase_array[i] = carg(X[i]) / (2 * M_PI);
      if (xkdb_array) xkdb_array[i] = 20 * log10(xkamp);
    }
    free(x);
  }
}

int total_interval_count = 0;

void format_0_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
                       int samples_per_sum, double begin_cutoff, double end_cutoff,
                       int interval_number,
                       int format_0_verbose,
                       int format_1_bins, int fundemental_db,
                       int call_format_1_analysis, int format_1_verbose,
                       int format_4_bins, double format_4_bin_size,
                       int call_format_4_analysis, int format_4_verbose)
{
  int interval_start = -1;
  int interval_end = -1;
  int interval_count = 0;
  if (last_frame == 0) last_frame = info->frames - 1;
  for (int frame = first_frame; frame <= last_frame; frame += samples_per_sum) {
    double sum = 0.0;
    int count =  ((last_frame - frame) > samples_per_sum) ? samples_per_sum : (last_frame - frame);
    for (int f = 0; f < count; f++) {
      double value = data[(frame+f)*info->channels+channel];
      sum += value*value;
    }
    double normalized_sum = sum/count;
    if (normalized_sum < ((interval_start == -1) ?  begin_cutoff : end_cutoff)) {
      if (interval_start > -1) {
        double interval_size = ((double)(1+interval_end-interval_start))/info->samplerate;
        if (format_0_verbose)
          printf("<<< end %d size=%.2f\n", interval_count, interval_size);
        else {
          printf("interval_number=%d duration=%.2f time=(%8.3f,%8.3f)\n", interval_count, interval_size,
                 ((double)interval_start)/info->samplerate, ((double)interval_end)/info->samplerate);
          fflush(stdout);
        }
        if (interval_size < 0.3) {
          interval_count --;
        } if (call_format_1_analysis && (interval_number < 0 || interval_number == interval_count)) {
          total_interval_count++;
          int first_frame = interval_start;
          int last_frame = frame;
          double discovered_note_number;
          format_1_analysis(info, data, channel, first_frame, last_frame, 
                            format_1_bins, fundemental_db, interval_count,
                            &discovered_note_number, format_1_verbose);
          fflush(stdout);
          if (call_format_4_analysis && !use_fft) {
            format_4_analysis(info, data, channel, first_frame, last_frame, format_4_bins, format_4_bin_size,
                              discovered_note_number, format_4_verbose);
            fflush(stdout);
            if (total_interval_count != lrint(discovered_note_number)) {
              printf("retry with interval_count as the note_number\n");
              format_4_analysis(info, data, channel, first_frame, last_frame, format_4_bins, format_4_bin_size,
                                (double)total_interval_count, format_4_verbose);
              fflush(stdout);
            }
          }
          printf("\n");
        }            
        interval_start = -1;
        interval_end = -1;
      }
    } else {
      if (interval_start == -1) {
        interval_count ++;
        interval_start = frame;
        if (format_0_verbose) printf(">>> begin %d\n", interval_count);
      }
      interval_end = frame;
    }
    if (format_0_verbose) {
      printf("%9.2f %8.4e\n", ((double)frame)/info->samplerate, normalized_sum);
    }
  }
}

#define D_INDEX_COUNT 10
                   
void show_change(FILE *out, double *note_number_array, double *frequency_array, double *xkdb_array,
                 int result_bins, int index, double xkdb_max)
{
  if (xkdb_array[index] < xkdb_array[index-1] || xkdb_array[index] < xkdb_array[index+1])
    return;
  if ((xkdb_max-xkdb_array[index]) > 30)
    return;
  double p = (1.0/2.0)*(xkdb_array[index-1]-xkdb_array[index+1])/(xkdb_array[index-1]+xkdb_array[index+1]-2*xkdb_array[index]);
  double xkdbp = xkdb_array[index] - (1.0/4.0)*(xkdb_array[index-1]-xkdb_array[index+1])*p;
  double frequencyp = frequency_array[index] + p*(frequency_array[index+1]-frequency_array[index-1])/2;
  double db[D_INDEX_COUNT];
  char low_greater_than_max[D_INDEX_COUNT] = {FALSE};
  char low_out_of_bounds[D_INDEX_COUNT] = {FALSE};
  int low_index[D_INDEX_COUNT] = {0};
  char high_greater_than_max[D_INDEX_COUNT] = {FALSE};
  char high_out_of_bounds[D_INDEX_COUNT] = {FALSE};
  int high_index[D_INDEX_COUNT] = {0};
  for (int d_index=0; d_index<D_INDEX_COUNT; d_index++) {
    double dbv = (d_index + 1) * 3.0;
    db[d_index] = dbv;
    int i;
    for (i=index-1; i>=0; i--) {
      if (xkdb_array[i] > xkdbp) {
        low_greater_than_max[d_index] = TRUE;
        break;
      }
      if (xkdb_array[i] < (xkdbp-dbv)) {
        break;
      }
    }
    if (i < 0) {
      low_out_of_bounds[d_index] = TRUE;
    }
    low_index[d_index] = i;
    for (i=index+1; i<result_bins; i++) {
      if (xkdb_array[i] > xkdbp) {
        high_greater_than_max[d_index] = TRUE;
        break;
      }
      if (xkdb_array[i] < (xkdbp-dbv)) {
        break;
      }
    }
    if (i < 0) {
      high_out_of_bounds[d_index] = TRUE;
    }
    high_index[d_index] = i;
  }
  if (use_formatted_show_change) {
    if (low_greater_than_max[6] || high_greater_than_max[6])
      return;
    fprintf(out, "%6.3f %12.4f %8.4f ", frequency_to_note_number(frequencyp), frequencyp, xkdbp-xkdb_max);
    for (int d_index=0; d_index<D_INDEX_COUNT; d_index++) {
      double dbv = db[d_index];
      if (low_greater_than_max[d_index] || high_greater_than_max[d_index])
        break;
      if (low_out_of_bounds[d_index] || high_out_of_bounds[d_index])
        break;
      if ((dbv+xkdb_max-xkdbp)>48)
        break;
      if (d_index > 0) fprintf(out, ", ");
      if (d_index % 4 == 0) fprintf(out, "\n       ");
      fprintf(out, "[%02.0f,%6.3f,%6.3f]", dbv, note_number_array[low_index[d_index]], note_number_array[high_index[d_index]]);
    }
    fprintf(out, "\n");
  } else {
    if (low_greater_than_max[6] || high_greater_than_max[6])
      return;
    int last_included_index = -1;
    for (int d_index=0; d_index<D_INDEX_COUNT; d_index++) {
      if (low_greater_than_max[d_index] || high_greater_than_max[d_index] ||
          low_out_of_bounds[d_index] || high_out_of_bounds[d_index])
        break;
      double dbv = db[d_index];
      if ((dbv+xkdb_max-xkdbp)>48)
        break;
      last_included_index = d_index;
    }
    int low = low_index[last_included_index];
    int high = high_index[last_included_index];
    fprintf(out, "%6.3f %8.3f\n", note_number_array[low-1], 0.0);
    for (int i=low; i<=high; i++) {
      fprintf(out, "%6.3f %8.3f\n", note_number_array[i], xkdb_array[i]);
    }
    fprintf(out, "%6.3f %8.3f\n", note_number_array[high+1], 0.0);
  }
}

void format_1_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, int fundemental_db, int interval_number,
		       double *note_number_ptr, int verbose)
{
  double *xkdb_array;
  double *note_number_array;
  double *frequency_array;
  int result_bins = bins;
  transformn(info, data, channel, first_frame, last_frame,
	     1.0 - 0.5, 1.0/bins, 88*bins, &result_bins,
	     &note_number_array, &frequency_array,
	     NULL, NULL, &xkdb_array);
  double xkdb_max = -1000;
  for (int index = 0; index<result_bins; index++) {
    double xkdb = xkdb_array[index];
    if (xkdb > xkdb_max) {
      xkdb_max = xkdb;
    }
    if (verbose) {
      char note_name[13];
      note_number_to_note_name(note_number_array[index], note_name, sizeof(note_name), TRUE);
      printf("%12.4f %8.4f %-10s %8.4f\n", frequency_array[index], note_number_array[index], note_name, xkdb);
    }
  }
  printf("xkdb_max=%8.3f\n", xkdb_max);
  char out_filename[64];
  snprintf(out_filename, sizeof(out_filename), "%s/interval_%d_format_1.data", directory, interval_number);
  FILE *out = fopen(out_filename, "w");
  if (out) {
    for (int index = 0; index<result_bins; index++) {
      fprintf(out, "%12.4f %8.4f\n", frequency_array[index], xkdb_array[index]-xkdb_max);
    }
    fclose(out);
  }
  if (use_fft) {
    snprintf(out_filename, sizeof(out_filename), "%s/interval_%d_format_1_peak.data", directory, interval_number);
    out = fopen(out_filename, "w");
    if (out) {
      if (!use_formatted_show_change) {
        fprintf(out, "%d\n", last_frame-first_frame+1);
      }
      for (int index = 1; index<(result_bins-1); index++) {
        show_change(out, note_number_array, frequency_array, xkdb_array, result_bins, index, xkdb_max);
      }
      fclose(out);
    }
  }
#define PEAK_LIMIT 30
  int peak_count = 0;
  int peak_index[PEAK_LIMIT];
  double xkdb_pmax = -1000;
  int index_pmax = -1;
  for (int index = 0; index<result_bins; index++) {
    double xkdb = xkdb_array[index];
    if (verbose) printf("  %d %8.4f %8.4f %d\n", index, xkdb, xkdb_pmax, index_pmax);
    if (xkdb > xkdb_pmax) {
      xkdb_pmax = xkdb;
      index_pmax = index;
      if (verbose) printf("+ %d %8.4f %8.4f %d\n", index, xkdb, xkdb_pmax, index_pmax);
    } else if ((xkdb_pmax - xkdb) >= 20) {
      if (verbose) printf("- %d %8.4f %8.4f\n", index_pmax, (xkdb_pmax - xkdb), (xkdb_max - xkdb_pmax));
      if (index_pmax>=0 && (xkdb_max - xkdb_pmax) <= 26 && peak_count < PEAK_LIMIT) {
        peak_index[peak_count++] = index_pmax;
        if (verbose) printf("peak %d\n", peak_count);
      }
      xkdb_pmax = -1000;
      index_pmax = -1;
      if (verbose) printf("- %d %8.4f %8.4f %d\n", index, xkdb, xkdb_pmax, index_pmax);
    }
  }
  double diff_freq_array[PEAK_LIMIT];
  int qp = -1;
  for (int p = 0; p < peak_count; p++) {
    double xkdb_rel = xkdb_array[peak_index[p]] - xkdb_max;
    double frequency = frequency_array[peak_index[p]];
    double note_number = note_number_array[peak_index[p]];
    double diff_freq = p ? (frequency_array[peak_index[p]] - frequency_array[peak_index[p-1]]) : 100000.0;
    if (p>0) diff_freq_array[p-1] = diff_freq;
    if ((p==2 || diff_freq_array[p-3] > 20) && diff_freq_array[p-2] > 20 && diff_freq_array[p-1] > 20 && qp == -1) qp = p-2;
    char note_name[13];
    note_number_to_note_name(note_number, note_name, sizeof(note_name), TRUE);
    if (p==0)
      printf("frequency=%9.3f,                     xkdb_rel=%8.4f, name=%s\n", frequency, xkdb_rel, note_name);
    else
      printf("frequency=%9.3f, diff_freq=%8.3f, xkdb_rel=%8.4f, name=%s\n", frequency, diff_freq, xkdb_rel, note_name);
  }
  if (note_number_ptr) {
    if (peak_count <= 8) {
      printf("** %8.2f\n", (xkdb_array[peak_index[1]] - xkdb_array[peak_index[0]]));
      if (peak_count < 2 || (xkdb_array[peak_index[1]] - xkdb_array[peak_index[0]]) <= 12) /* 18 */
        *note_number_ptr = note_number_array[peak_index[0]];
      else
        *note_number_ptr = note_number_array[peak_index[1]];
    } else {
      double note_number_qp = frequency_to_note_number(diff_freq_array[qp]);
      *note_number_ptr = note_number_qp - 0.05 * (frequency_array[peak_index[qp+1]]/diff_freq_array[qp]);
    }
  }
  free(xkdb_array);
  free(note_number_array);
  free(frequency_array);
}

void format_2_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double)
{
  int note_number = lrint(note_number_double);
  double note = 440 * exp(log(2.0) * (note_number - 49) / 12);
  for (int harmonic=1; harmonic<=8; harmonic++) {
    printf("harmonic=%d\n", harmonic);
    for (int b = 0; b<bins; b++) {
      double bin = bin_size*((double)b/bins-0.5);
      double freq = note * harmonic * exp(log(2.0) * bin/12);
      double xkamp, xkphase, xkdb;
      transform1(info, data, channel, first_frame, last_frame,
		 freq, &xkamp, &xkphase, &xkdb);
      printf("bin=%7.2f, freq=%8.4e, 20*log10(amp)=%8.5f, amp=%8.5e, phase=%8.5f\n", bin*100, freq, xkdb, xkamp, xkphase);
    }
    printf("\n");
  }
}

void format_3_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double)
{
  int note_number = lrint(note_number_double);
  double note = 440 * exp(log(2.0) * (note_number - 49) / 12);
  for (int b = 0; b<bins; b++) {
    double bin = bin_size*((double)b/bins-0.5);
    double fundamental = note * exp(log(2.0) * bin / 12);
    printf("%8.3f, %9.3f", bin*100, fundamental);
    for (int harmonic=1; harmonic<=8; harmonic++) {
      double freq = harmonic * fundamental;
      double xkdb;
      transform1(info, data, channel, first_frame, last_frame,
		 freq, NULL, NULL, &xkdb);
      printf(", %8.5f", xkdb);
    }
    printf("\n");
  }
}

void format_4_analysis(SF_INFO *info, double *data, int channel, int first_frame, int last_frame,
		       int bins, double bin_size,
		       double note_number_double, int verbose)
{
  int note_number = lrint(note_number_double);
  double note_frequency = note_number_to_frequency(note_number);
  printf("note_number=%d, frequency=%9.2f\n", note_number, note_number_to_frequency((double)note_number));
  double *xkdb_array;
  double *note_number_array;
  int max_harmonics = 1 + lrint(((note_number < 13) ? 1000 :
                                 (note_number < 25) ? 2000 :
                                 (note_number < 37) ? 4000 : 8000)/note_frequency);
  for (int harmonic=1; harmonic<=max_harmonics; harmonic++) {
    double harmonic_note_number = note_number + 12 * log((double)harmonic) / log(2);
    int result_bins = bins;
    transformn(info, data, channel, first_frame, last_frame,
               harmonic_note_number - bin_size/2, bin_size/bins, bins, &result_bins,
               &note_number_array, NULL,
               NULL, NULL, &xkdb_array);
    double xkdb_max = -10000;
    for (int b = 0; b<bins; b++) {
      if (xkdb_max < xkdb_array[b]) xkdb_max = xkdb_array[b];
    }
    if (xkdb_max < -90) continue;
    double xkdb_m20 = xkdb_max - 20; int xkdb_m20_first = -1; int xkdb_m20_last = -1;
    double xkdb_m12 = xkdb_max - 12; int xkdb_m12_first = -1; int xkdb_m12_last = -1;
    double xkdb_m06 = xkdb_max - 06; int xkdb_m06_first = -1; int xkdb_m06_last = -1;
    double xkdb_m00 = xkdb_max - 00; int xkdb_m00_first = -1; int xkdb_m00_last = -1;
    for (int b = 0; b<bins; b++) {
      if (xkdb_m20 <= xkdb_array[b]) {
	if (xkdb_m20_first < 0) xkdb_m20_first = b;
	xkdb_m20_last = b;
      }
      if (xkdb_m12 <= xkdb_array[b]) {
	if (xkdb_m12_first < 0) xkdb_m12_first = b;
	xkdb_m12_last = b;
      }
      if (xkdb_m06 <= xkdb_array[b]) {
	if (xkdb_m06_first < 0) xkdb_m06_first = b;
	xkdb_m06_last = b;
      }
      if (xkdb_m00 <= xkdb_array[b]) {
	if (xkdb_m00_first < 0) xkdb_m00_first = b;
	xkdb_m00_last = b;
      }
    }
    double note_number_max = note_number_array[(xkdb_m00_first+xkdb_m00_last)/2];
    double frequency_max = note_number_to_frequency(note_number_max);
    char note_name_max[13];
    note_number_to_note_name(note_number_max, note_name_max, sizeof(note_name_max), TRUE);
#define BIN(b) bin_size*((double)b/bins-0.5)
    char xkdb_m00_last_string[32];
    if (xkdb_m00_first == xkdb_m00_last)
      strcpy(xkdb_m00_last_string, "");
    else
      snprintf(xkdb_m00_last_string, sizeof(xkdb_m00_last_string), ",%7.2f", 100*BIN(xkdb_m00_last));
    printf("[%2d %-10s %10.4f %8.4f 0:(%7.2f%s), 6:(%7.2f,%7.2f), 12:(%7.2f,%7.2f), 20:(%7.2f,%7.2f)]",
	   harmonic, note_name_max, frequency_max, xkdb_max,
	   100*BIN(xkdb_m00_first), xkdb_m00_last_string,
	   100*BIN(xkdb_m06_first), 100*BIN(xkdb_m06_last),
	   100*BIN(xkdb_m12_first), 100*BIN(xkdb_m12_last),
	   100*BIN(xkdb_m20_first), 100*BIN(xkdb_m20_last));
    printf("\n");
    fflush(stdout);
  }
  free(xkdb_array);
  free(note_number_array);
}


int main(int argc, char **argv)
{
  char *filename = NULL;
  int channel = 0;
  int format = 0;
  int first_frame = 0;
  int last_frame= 0;
  double first_time = 0.0;
  double last_time = 0.0;
  int interval_number = -1;
  char *note_name = NULL;
  double note_number = 0.0;
  int fundemental_db = -20;
  int format_1_bins = 100;
  int bins = 600;
  double bin_size = 1.5;
  int verbose = TRUE;
  double seconds_per_sum = 0.05;
  double begin_cutoff = 1e-3;
  double end_cutoff = 1e-4;
  int format_0_verbose = FALSE;
  int call_format_1_analysis = TRUE;
  int format_1_verbose = FALSE;
  int call_format_4_analysis = TRUE;
  int format_4_verbose = FALSE;

  for (int p = 1; p < argc; p++) {
    if (0 == strcmp(argv[p], "-f") || 0 == strcmp(argv[p], "--file")) {
      if (++p < argc) filename = argv[p];
    } else if (0 == strcmp(argv[p], "-c") || 0 == strcmp(argv[p], "--channel")) {
      if (++p < argc) channel = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-fm") || 0 == strcmp(argv[p], "--format")) {
      if (++p < argc) format = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-ff") || 0 == strcmp(argv[p], "--firstframe")) {
      if (++p < argc) first_frame = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-lf") || 0 == strcmp(argv[p], "--lastframe")) {
      if (++p < argc) last_frame = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-ft") || 0 == strcmp(argv[p], "--firsttime")) {
      if (++p < argc) first_time = strtod(argv[p], NULL);
    } else if (0 == strcmp(argv[p], "-lt") || 0 == strcmp(argv[p], "--lasttime")) {
      if (++p < argc) last_time = strtod(argv[p], NULL);
    } else if (0 == strcmp(argv[p], "-i") || 0 == strcmp(argv[p], "--intervalnumber")) {
      if (++p < argc) interval_number = strtol(argv[p], NULL, 10);
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
    } else if (0 == strcmp(argv[p], "--format_1_bins")) {
      if (++p < argc) format_1_bins = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-b") || 0 == strcmp(argv[p], "--bins")) {
      if (++p < argc) bins = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-f1b") || 0 == strcmp(argv[p], "--format_1_bins")) {
      if (++p < argc) format_1_bins = strtol(argv[p], NULL, 10);
    } else if (0 == strcmp(argv[p], "-v") || 0 == strcmp(argv[p], "--verbose")) {
      if (++p < argc) verbose = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-f0v") || 0 == strcmp(argv[p], "--format_0_verbose")) {
      if (++p < argc) format_0_verbose = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-f1v") || 0 == strcmp(argv[p], "--format_1_verbose")) {
      if (++p < argc) format_1_verbose = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-f4v") || 0 == strcmp(argv[p], "--format_4_verbose")) {
      if (++p < argc) format_4_verbose = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-cf1") || 0 == strcmp(argv[p], "--call_format_1")) {
      if (++p < argc) call_format_1_analysis = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-cf4") || 0 == strcmp(argv[p], "--call_format_4")) {
      if (++p < argc) call_format_4_analysis = 0 == strcmp(argv[p], "true");
    } else if (0 == strcmp(argv[p], "-t1") || 0 == strcmp(argv[p], "--use_transform1")) {
      use_transform1 = TRUE;
    } else if (0 == strcmp(argv[p], "-fft") || 0 == strcmp(argv[p], "--use_fft")) {
      use_transform1 = FALSE;
      use_fft = TRUE;
    } else if (0 == strcmp(argv[p], "-cz") || 0 == strcmp(argv[p], "--use_chirpz")) {
      use_transform1 = FALSE;
      use_fft = FALSE;
      use_chirp_z = TRUE;
    } else {
      printf("Option %s was seen, but it is not supported\n", argv[p]);
    }
  }
  if (filename == NULL) {
    printf("Please specify a filename, using -f\n");
    exit(0);
  }

  directory = use_transform1 ? "data_transform1" : use_fft ? "data_fft" : "data_chirpz";
  double note = 440 * exp(log(2.0) * note_number / 12);
  if (format > 1) {
    printf("note_name=%s, note_number %8.3f, freq=%9.3f\n", note_name, note_number, note);
  }
  
  SF_INFO info_s, *info = &info_s;
  memset(info, 0, sizeof(*info));
  SNDFILE *in = sf_open(filename, SFM_READ, info);
  printf("format=%08X, frames=%d, samplerate=%d, channels=%d, sections=%d\n",
	 info->format, (int)info->frames, info->samplerate, info->channels, info->sections);
  int data_count = info->frames * info->channels;
  double *data = (double *)calloc(data_count, sizeof(double));
  int read_count = sf_read_double(in, data, data_count);
  if (data_count != read_count) {
    printf("Unexpected result from sf_read_double\n");
  }
  sf_close(in);

  if (first_time != 0.0) first_frame = lrint(first_time * info->samplerate);
  if (last_time != 0.0) last_frame = lrint(last_time * info->samplerate);
  int samples_per_sum = (int)(seconds_per_sum * info->samplerate);
  
  clock_t start_clock = clock();
  if (format == 0) {
    format_0_analysis(info, data, channel, first_frame, last_frame,
                      samples_per_sum, begin_cutoff, end_cutoff,
                      interval_number,
                      format_0_verbose,
                      format_1_bins, fundemental_db,
                      call_format_1_analysis, format_1_verbose,
                      bins, bin_size,
                      call_format_4_analysis, format_4_verbose);
  } else if (format == 1) {
    double discovered_note_number;
    format_1_analysis(info, data, channel, first_frame, last_frame, format_1_bins,
                      fundemental_db, 0,
		      &discovered_note_number, verbose);
  } else if (format == 2) {
    format_2_analysis(info, data, channel, first_frame, last_frame, bins, bin_size,
		      note_number);
  } else if (format == 3) {
    format_3_analysis(info, data, channel, first_frame, last_frame, bins, bin_size,
		      note_number);
  } else if (format == 4) {
    format_4_analysis(info, data, channel, first_frame, last_frame, bins, bin_size,
		      note_number, verbose);
  }
  clock_t end_clock = clock();
  double cpu_time = ((double)(end_clock-start_clock))/CLOCKS_PER_SEC;
  printf("\n\ntotal cpu time=%9.2f, intervals=%d\n", cpu_time, total_interval_count);
  double transform1_cpu_time = ((double)transform1_clock)/CLOCKS_PER_SEC;
  printf("transform1 cpu_time=%9.2f, calls=%lld, frames=%8.3e\n", transform1_cpu_time, transform1_calls, (double)transform1_frames);
}
