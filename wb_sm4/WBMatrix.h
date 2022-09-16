#ifndef _HWBMATRIX_H_
#define _HWBMATRIX_H_

#include "structure.h"
#include "random.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>



void SetRandSeed(unsigned int seed);//Set random seed
int xorU8(uint8_t n);
int xorU16(uint16_t n);
int xorU32(uint32_t n);
//void V32toV4(uint32_t n, uint8_t* b[8]);
//void V4toV32(uint32_t* n, uint8_t b[8]);
void initV8(V8* Vec);
void initV32(V32* Vec);
void MatMulVecM32(M32 Mat, V32 Vec, V32* ans);
void MatMulVecM8(M8 Mat, V8 Vec, V8* ans);
uint32_t MatMulNumM32(M32 Mat, uint32_t n);   //32*32矩阵与32维向量相乘
uint32_t VecAddVecV32(uint32_t Vec1, uint32_t Vec2);   //32bit向量的异或
void identityM8(M8* Mat);
void identityM32(M32* Mat);
void randV8(V8* Vec);
void randV32(V32* Vec);
void randM8(M8* Mat);
void randM32(M32* Mat);
void copyM8(M8 Mat1, M8* Mat2);
void copyM32(M32 Mat1, M32* Mat2);
int isequalM32(M32 Mat1, M32 Mat2);
int isinvertM32(M32 Mat);
uint32_t affineU32(Aff32 aff, uint32_t arr);
void genMatpairM8(M8* Mat, M8* Mat_inv);
void genMatpairM32(M32* Mat, M32* Mat_inv);
void genaffinepairM8(Aff8* aff, Aff8* aff_inv);
void genaffinepairM32(Aff32* aff, Aff32* aff_inv);
void affinecomM8to32(Aff8 aff1, Aff8 aff2, Aff8 aff3, Aff8 aff4, Aff32* aff);
void MatrixcomM8to32(M8 m1, M8 m2, M8 m3, M8 m4, M32* mat);
void VectorcomV8to32(V8 v1, V8 v2, V8 v3, V8 v4, V32* vec);

#endif