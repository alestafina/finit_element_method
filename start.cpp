#include "fem.h"

int main() {
   CGM matrix;
   matrix.read_data();
   matrix.glob_matrix();
   matrix.read_boundary();
   //matrix.thrid_cond();
   //matrix.second_cond();
   matrix.first_cond();
   matrix.CGM_precond_ILU();
   matrix.print_result();

   return 0;
}
