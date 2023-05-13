#pragma once

#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

struct matrix {
   vector<int> ig, jg;
   vector<double> al, au, di;
};

class FEM {
protected:
   matrix A;                        // �������� �������
   vector<vector<double>> loc_M;    // ��������� ������� ����
   vector<vector<double>> loc_G;    // ��������� ������� ���������
   vector<double> loc_f;            // ��������� ������ ������ �����
   vector<double> F;                // ������ ������ �����
   vector<double> q;                // ������ ����� (�����)
   vector<double> q_real;           // ��������� �������� �����

   int hx, hy;
   vector<double> x_nodes, y_nodes; // ����� �� x � y
   vector<int> nodes;               // ���������� ������ ����� 
   vector<vector<int>> W;           // ������ "������" ������ ������� ��

   vector<int> first_bc, second_bc, thrid_bc; // ����, � ������� ��������� ��

   int N, M, nx, ny;                // ���-�� �������� ���������, ���-�� �����
   double betta, lambda, gamma;
public:
   FEM();
   ~FEM();

   void read_data();                // ���� ������
   void making_grid();              // ����������� �����
   void build_portrait();           // ����������� �������� ������� �� �����
   void glob_matrix();              // ������ ���������� �������
   void add_elem(int i, int j, double elem); // ���������� �������� � ���������� �������
   
   void read_boundary();            // ���� ������� �������
   void thrid_cond();               // ���� �������� �������� �������
   void second_cond();              // ���� ������� �������� �������
   void first_cond();               // ���� ������� �������� �������

   void print_result();
};

class Functions {
public:
   double func(double x, double y);
   double u_g(double x, double y);
   double u_betta(double x, double y);
   double theta(double x, double y);
   double u_real(double x, double y);
};

class CGM : public FEM {
private:
   vector<double> r, z, buf1, buf2;
   vector<double> al_LU, au_LU, di_LU;
   int maxiter, k;
   double eps;
   double residual, norm;
public: 
   CGM();
   ~CGM();

   void ILU_precond(); // �������� LU ������������
   
   void set_vector();
   vector<double> matr_vec_mult(vector<double> &x, bool flag);
   double dot_product(vector<double> &a, vector<double> &b);

   vector<double> direct(vector<double> &L, vector<double> &b);
   vector<double> direct(vector<double> &L, vector<double> &D, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &b);
   vector<double> reverse(vector<double> &U, vector<double> &D, vector<double> &b);
   
   void CGM_precond_ILU();
};
