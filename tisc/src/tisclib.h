#ifndef TISCLIB_H
#define TISCLIB_H

struct GRIDNODE;
struct BLOCK;

int Allocate_Memory_for_external_use();
int Allocate_Memory();
float calculate_topo(ModelConfig *cfg, ModelContext *ctx, float **topo_new);
int Delete_Block(int i_Block);
int insert_new_Block(int num_new_Block);
int gradual_Block(ModelConfig *cfg, ModelContext *ctx);
int defineLESalmostdiagonalmatrix (ModelConfig *cfg, double **A, double *b, float **q, float **Dq, float **w, bool doing_visco);
int defineLESmatrix_for_mathlib (ModelConfig *cfg, float *A, int *IA, int *JA, float *b, float **q, float **w, int *nonzeroes, bool doing_visco);
int Matrix_Info(float **xarxa, int nx, int ny, float dx, float dy, float *max, float *min, float *maxcurvx, float *maxcurvy, float *maxcurvxy);
int Perfil_info(float *perfil, int n, float *max, float *min);
int Repare_Blocks(ModelConfig *cfg, ModelContext *ctx);
int solveLESalmostdiagonal (ModelConfig *cfg, double **A, double *b, float **x);
float Sort_Matrix (float **matrix, struct GRIDNODE *orden, int Nx, int Ny);
float ReSort_Matrix (float **matrix, struct GRIDNODE *orden, int Nx, int Ny);
float calculate_sea_level(float current_time);
int calculate_water_load(ModelConfig *cfg, ModelContext *ctx);

#endif /* TISCLIB_H */