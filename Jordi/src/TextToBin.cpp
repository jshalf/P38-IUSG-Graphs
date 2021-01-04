#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int main(int argc, char *argv[])
{
   int arg_iter = 0;
   char in_filename[1000], out_filename[1000];
   int file_flag = 0, dir_flag = 0;
   int row, col;
   double elem;
   size_t size;

   while (arg_iter < argc){
      if (strcmp(argv[arg_iter], "-file") == 0){
         arg_iter++;
         strcpy(in_filename, argv[arg_iter]);
         file_flag = 1;
      }
      else if (strcmp(argv[arg_iter], "-dir") == 0){
         dir_flag = 1;
      }
      arg_iter++;
   }

   if (file_flag){
      sprintf(out_filename, "%s.bin", in_filename);
      FILE *out_file = fopen(out_filename, "wb");
      FILE *in_file = fopen(in_filename, "r");
      while(fscanf(in_file, "%d %d %lg", &row, &col, &elem) == 3){
         fwrite(&row, sizeof(int), 1, out_file);
         fwrite(&col, sizeof(int), 1, out_file);
         fwrite(&elem, sizeof(double), 1, out_file);
      }
      fclose(out_file);
      fclose(in_file);
   } 
   return 0;
}
