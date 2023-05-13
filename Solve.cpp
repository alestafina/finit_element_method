#include "fem.h"

CGM::CGM() {
   r = vector<double>();
   z = vector<double>();
   q = vector<double>();
   buf1 = vector<double>();
   buf2 = vector<double>();
   al_LU = vector<double>();
   au_LU = vector<double>();
   di_LU = vector<double>();
   maxiter = 1000;
   k = 0;
   eps = 1e-20;
   residual = 10.0;
   norm = 0;
}

CGM::~CGM() {
   r.~vector();
   z.~vector();
   q.~vector();
   buf1.~vector();
   buf2.~vector();
   al_LU.~vector();
   au_LU.~vector();
   di_LU.~vector();
}

void CGM::ILU_precond() {
   al_LU = A.al;
   au_LU = A.au;
   di_LU = A.di;
   double sumL = 0.0, sumU = 0.0, sumD = 0.0;
   for (int k = 0; k < M; k++) {
      sumD = 0.0;
      int i0 = A.ig[k], i1 = A.ig[k + 1];
      int i = i0;
      for (; i0 < i1; i0++) {
         sumL = 0.0;
         sumU = 0.0;
         int j0 = i, j1 = i0;
         for (; j0 < j1; j0++) {
            int t0 = A.ig[A.jg[i0]], t1 = A.ig[A.jg[i0] + 1];
            for (; t0 < t1; t0++) {
               if (A.jg[j0] == A.jg[t0]) {
                  sumL += al_LU[j0] * au_LU[t0];
                  sumU += au_LU[j0] * al_LU[t0];
               }
            }
         }
         al_LU[i0] -= sumL;
         au_LU[i0] -= sumU;
         au_LU[i0] /= di_LU[A.jg[i0]];
         sumD += al_LU[i0] * au_LU[i0];
      }
      di_LU[k] -= sumD;
   }
}

void CGM::set_vector() {
   q.resize(M);
}

vector<double> CGM::matr_vec_mult(vector<double> &x, bool flag) {
   vector<double> result(x.size()), L = A.al, U = A.au;
   if (flag) swap(L, U);
   for (int i = 0; i < x.size(); ++i) {
      result[i] = A.di[i] * x[i];
      for (int j = A.ig[i]; j < A.ig[i + 1]; j++) {
         result[i] += L[j] * x[A.jg[j]];
         result[A.jg[j]] += U[j] * x[i];
      }
   }
   return result;
}

double CGM::dot_product(vector<double> &a, vector<double> &b) {
   double result = 0.0;
   for (int i = 0; i < N; ++i)
      result += a[i] * b[i];
   return result;
}

vector<double> CGM::direct(vector<double> &L, vector<double> &b) { // прямой U(T)x = b
   vector<double> result(M);
   for (int i = 0; i < M; i++) {
      double sum = 0.0;
      for (int k0 = A.ig[i]; k0 < A.ig[i + 1]; k0++) {
         sum += result[A.jg[k0]] * L[k0];
      }
      result[i] = b[i] - sum;
   }
   return result;
}

vector<double> CGM::direct(vector<double> &L, vector<double> &D, vector<double> &b) { // пря-мой Lx = b
   vector<double> result(M);
   for (int i = 0; i < M; i++) {
         double sum = 0.0;
      for (int k0 = A.ig[i]; k0 < A.ig[i + 1]; k0++) {
         sum += result[A.jg[k0]] * L[k0];
      }
      result[i] = (b[i] - sum) / D[i];
   }
   return result;
}

vector<double> CGM::reverse(vector<double> &U, vector<double> &b) { // обратный Ux = b
   vector<double> result = b;
   for (int i = M - 1; i >= 0; i--) {
      int k0 = A.ig[i], k1 = A.ig[i + 1];
      int j;
      for (; k0 < k1; k0++) {
         j = A.jg[k0];
         result[j] -= result[i] * U[k0];
      }
   }
   return result;
}

vector<double> CGM::reverse(vector<double> &U, vector<double> &D, vector<double> &b) { // об-ратный L(T)x = b
   vector<double> result = b;
   for (int i = M - 1; i >= 0; i--) {
      int k0 = A.ig[i], k1 = A.ig[i + 1];
      int j;
      result[i] /= D[i];
      for (; k0 < k1; k0++) {
         j = A.jg[k0];
         result[j] -= result[i] * U[k0];
      }
   }
   return result;
}

void CGM::CGM_precond_ILU() {
   double a, b, r_next, r_prev;

   set_vector();
   ILU_precond();

   norm = dot_product(F, F);
   buf1 = matr_vec_mult(q, 0);
   
   r.resize(M);
   
   for (int i = 0; i < M; i++) { // r0 = f - Ax0
      r[i] = F[i] - buf1[i];
   }
   r = direct(al_LU, di_LU, r); // r0 = L(-1)(f - Ax0)
   r = reverse(al_LU, di_LU, r); // r0 = L(-T)L(-1)(f - Ax0)
   r = matr_vec_mult(r, 1); // r0 = A(T)L(-T)L(-1)(f - Ax0)
   r = direct(au_LU, r); // r0 = U(-T)A(T)L(-T)L(-1)(f - Ax0)
   z = r;
   
   buf1.clear();
   buf1.resize(M);
   for (int i = 0; i < M; i++) { // x = Ux
      buf1[i] += q[i];
      for (int j = A.ig[i]; j < A.ig[i + 1]; j++)
         buf1[A.jg[j]] += au_LU[j] * q[i];
   }
   q = buf1;

   while (k < maxiter && residual > eps) {
      r_prev = dot_product(r, r); // r(k-1)
      buf2 = reverse(au_LU, z); // U(-1)z
      buf2 = matr_vec_mult(buf2, 0); // AU(-1)z
      buf2 = direct(al_LU, di_LU, buf2); // L(-1)AU(-1)z
      buf2 = reverse(al_LU, di_LU, buf2); // L(-T)L(-1)AU(-1)z
      buf2 = matr_vec_mult(buf2, 1); // A(T)L(-T)L(-1)AU(-1)z
      buf2 = direct(au_LU, buf2); // U(-T)A(T)L(-T)L(-1)AU(-1)z
      a = r_prev / dot_product(buf2, z);
      for (int i = 0; i < M; ++i) {
         q[i] += a * z[i];
         r[i] -= a * buf2[i];
      }
      r_next = dot_product(r, r);
      b = r_next / r_prev;
      for (int i = 0; i < M; ++i) {
         z[i] = r[i] + b * z[i];
      }
      residual = sqrt(r_next / norm);
      ++k;
   }
   q = reverse(au_LU, q); // x = U(-1)x
}
