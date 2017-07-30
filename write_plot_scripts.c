/*
gcc -g -O2 -std=gnu99 -Wall -o write_plot_scripts write_plot_scripts.c 
*/

#include <stdio.h>

int main(int argc, char **argv)
{
  printf("set term png size 16000, 8000\n");
  printf("set logscale x\n");
  printf("set xrange [25:5000]\n");
  printf("set yrange [-90:0]\n");
  printf("set grid xtics\n");

  for (int i = 1; i <= 88; i++) {
    printf("set output \"data/interval_%d_format_1.png\"\n", i);
    printf("plot \"data/interval_%d_format_1.data\" using 1:2 with lines\n", i);
  }
}
