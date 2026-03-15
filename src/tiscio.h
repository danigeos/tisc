#ifndef TISCIO_H
#define TISCIO_H

#include <stdio.h>

struct BLOCK;
struct CS2D;

int reformat_file_thin_sheet_BC(ModelConfig *cfg, ModelContext *ctx, char *tmpTSBCfilename);
int read_file_horiz_record_time();
int read_file_initial_deflection (float **w);
int read_file_initial_topo(float **h);
int read_file_initial_rivers();
int read_file_parameters (int show, int reformat);
int match_parameter (char *str1, char *str2, int show, int replace, char *line);
int read_file_rain(float **precipitation_file);
int read_file_resume(char *filename);
int read_file_insolation();
int read_file_sea_level();
int read_file_node_defs(ModelConfig *cfg, ModelContext *ctx, float dt_st);
int read_file_Te();
int read_file_output_Blocks ();
int read_file_2D_CS (struct BLOCK *Blocks, struct CS2D *CrossSection, int Nx2D);
int Calculate_2D_Cross_Section (ModelConfig *cfg, ModelContext *ctx, struct BLOCK *Blocks, struct CS2D *CrossSection, int Nx2D);
int write_file_cross_section (ModelConfig *cfg, ModelContext *ctx);
int write_file_Blocks (ModelConfig *cfg, ModelContext *ctx);
int write_file_grainsize (ModelConfig *cfg, ModelContext *ctx);
int write_file_thicksalt (ModelConfig *cfg, ModelContext *ctx);
int write_file_drainage (ModelConfig *cfg, ModelContext *ctx);
int write_file_river_basins (ModelConfig *cfg, ModelContext *ctx);
int write_file_lakes (ModelConfig *cfg, ModelContext *ctx);
int write_file_ice (ModelConfig *cfg, ModelContext *ctx);
int write_file_resume(ModelConfig *cfg, ModelContext *ctx);
int write_file_surftransp (ModelConfig *cfg, ModelContext *ctx);
int write_file_time (ModelConfig *cfg, ModelContext *ctx, float **xarxa1, float **xarxa2);
int write_file_deflection (ModelConfig *cfg, ModelContext *ctx);
int write_file_velocity_field (ModelConfig *cfg, ModelContext *ctx);

#endif /* TISCIO_H */