#include "wbsm4.h"

void printstate(unsigned char* in)
{
    int i;
    for (i = 0; i < 16; i++)
    {
        printf("%.2X", in[i]);
    }
    printf("\n");
}
void randomOutIn(unsigned char Out[32][4][8][16], unsigned char In[32][4][8][16]) {//只第一个可为空

    //初始化
    srand((unsigned int)time(NULL));		//时间播种
    for (int r = 0; r < 32; r++)					//前9轮
        for (int i = 0; i < 4; i++) {			//每一轮有4次列混合
            for (int j = 0; j < 8; j++) {		//每一次列混合，能用到的Out总数
                for (int k = 0; k < 16; k++) {	//大小为16的一维表，作用是4bit到4bit的随机代换
                    Out[r][i][j][k] = k;

                }
            }
        }

    //随机置换生成Out 
    for (int r = 0; r < 32; r++)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                for (int k = 0; k < 16; k++) {
                    int randNum = rand() % 16;
                	
                    int t = Out[r][i][j][k];
                    Out[r][i][j][k] = Out[r][i][j][randNum];
                    Out[r][i][j][randNum] = t;
                }
            }
        }

    //生成In 
    for (int r = 0; r < 32; r++)
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 8; j++) {
                for (int k = 0; k < 16; k++) {
                    In[r][i][j][Out[r][i][j][k]] = k;
                }
            }
        }
    return;
}

void split32to8_mat(Aff32 MB, Aff8* MB_ij[4][4]) { 
     for (int r = 0; r < 32; r++) {
          uint32_t m = MB.Mat.M[r]; 
          MB_ij[r / 8][0]->Mat.M[r % 8] = (m << 24) & 0xff;
          MB_ij[r / 8][1]->Mat.M[r % 8] = (m << 16) & 0xff;
          MB_ij[r / 8][2]->Mat.M[r % 8] = (m << 8) & 0xff;
          MB_ij[r / 8][3]->Mat.M[r % 8] = (m << 0) & 0xff;    
     }
     for (int i = 0; i < 4; i++) {
         for (int j = 0; j < 3; j++) {
             MB_ij[i][j]->Vec.V = 0;
         }
     }
     for (int i = 0; i < 4; i++) {
         MB_ij[i][3]->Vec.V = MB.Vec.V << (8 * (3 - i)) & 0xff;
     }
    
}


void wbsm4_gen_part1Table() {

}

void wbsm4_gen_part2Table(uint8_t* key)
 {
    sm4_context ctx;
    sm4_setkey_enc(&ctx, key);   //密钥扩展
    InitRandom(((unsigned int)time(NULL)));
    //生成32轮，每轮4个的  随机布尔掩码
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            Mask[i][j] = cus_random();
        }
    }
    for (int i = 0; i < 32; i++)
    {
        //affine E 
        for (int j = 0; j < 4; j++)
        {
            genaffinepairM8(&Eij[i][j], &Eij_inv[i][j]);
        }

        // combine 4 E8 to 1 E32
        affinecomM8to32(Eij[i][0], Eij[i][1], Eij[i][2], Eij[i][3], &Ei[i]);
    }


    //生成32轮，每轮四个的  仿射线性编码矩阵和常数以及其逆矩阵
    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            genaffinepairM32(&MB[i][j], &MB_inv[i][j]);
        }
    }
    //生成out以及in的非线性编码
    randomOutIn(Out_part2, In_part2);
    for (int i = 0; i < 4; i++)                         //第一轮进入part2不带编码，所以输入是x输出也是x
        for (int x = 0; x < 256; x++)
            Table_before_part2[0][i][x] = x;

    for (int r = 1; r < 32; r++) {                                //进入s盒之前，去掉part1的非线性编码
        for(int i = 0; i < 4; i++){
            for (int x = 0; x < 256; x++) {
                uint8_t y = x;
                Table_before_part2[r][i][x] =
                    (In_part1[r][2 * i + 0][y >> 4 & 0xf] << 4) |
                    (In_part1[r][2 * i + 1][y >> 0 & 0xf] << 0);

            }
        }
    }

    for (int r = 1; r < 32; r++) {                             //进入s盒之前，去掉上一轮的线性编码
        for (int i = 0; i < 4; i++) {
            for (int x = 0; x < 256; x++) {
                Table_before_part2[r][i][x] = affineU8(Eij_inv[r][i], Table_before_part2[r][i][x]);
            
            }
        }
    }

    for (int r = 1; r < 32; r++) {                         //进入s盒之前，去掉上一轮的掩码
        for (int i = 0; i < 4; i++) {
            for (int x = 0; x < 256; x++) {
                Table_before_part2[r][i][x] = Table_before_part2[r][i][x] ^ VecAddVecV8(Mask[r][0], Mask[r][1], Mask[r][2], Mask[r][3]);
            }
        }
    }                     

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
                uint8_t temp_u8 = SBOX[Table_before_part2[i][j][x] ^ (ctx.sk[i] >> (24 - j * 8)) && 0xFF];   //异或rk，sbox
                uint32_t temp_u32 = temp_u8 << (24 - j * 8);                      //8bit扩展成32bit
                Table_SL[i][j][x] = MatMulNumM32(L_matrix, temp_u32);             //异或rk，sbox，L移位
            }
        }
    }

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
               
                uint32_t temp_u32 = VecAddVecV32(Table_SL[i][j][x], Mask[i][j]);      //异或rk，sbox，L移位,+mask
                temp_u32 = MatMulNumM32(MB[i][j].Mat, temp_u32);    //异或rk，sbox，L移位,+mask,+线性矩阵MB.Mat
                Table_part2[i][j][x] = VecAddVecV32(temp_u32, MB[i][j].Vec.V);    //异或rk，sbox，L移位,+mask,+线性矩阵MB.Mat,+MB.Vec

            }
        }
    }
    for (int r = 0; r < 32; r++) {                                ////异或rk，sbox，L移位,+mask,+线性矩阵MB.Mat,+MB.Vec,+非线性编码Out
        for (int i = 0; i < 4; i++) {
            for (int x = 0; x < 256; x++) {
                uint8_t y = x;
                uint32_t a = Table_part2[r][i][y];
                Table_part2[r][i][x] = 
                    (Out_part2[r][i][0][a >> 28 & 0xf] << 28) |
                    (Out_part2[r][i][1][a >> 24 & 0xf] << 24) |
                    (Out_part2[r][i][2][a >> 20 & 0xf] << 20) |
                    (Out_part2[r][i][3][a >> 16 & 0xf] << 16) |
                    (Out_part2[r][i][4][a >> 12 & 0xf] << 12) |
                    (Out_part2[r][i][5][a >> 8 & 0xf] << 8)   |
                    (Out_part2[r][i][6][a >> 4 & 0xf] << 4)   |
                    (Out_part2[r][i][7][a >> 0 & 0xf] << 0);

            }
        }
    }
    //做一个8bit输入，32bit输出的将线性编码MB换成8*8线性编码E的表，在异或表之前完成
    //
    

    for (int r = 0; r < 32; r++) {
        for (int i = 0; i < 4; i++) {
            Aff8 MB_ij[4][4];
            V8 temp0[4];
            V8 temp1[4];
            split32to8_mat(MB[r][i], &MB_ij);
            for (int j = 0; j < 4; j++) {
                for (int x = 0; x < 16; x++) {
                    for (int y = 0; y < 16; y++) {
                        V8 temp_u8;
                        temp_u8.V = (In_part2[r][i][j * 2][x] << 4) |
                                    (In_part2[r][i][j * 2 + 1][y]);
                        temp_u8.V ^= MB[r][i].Vec.V << (8 * (3 - j)) & 0xff;
                        MatMulVecM8(MB_ij[0][j].Mat, temp_u8, &temp0[0]);
                        MatMulVecM8(MB_ij[1][j].Mat, temp_u8, &temp0[1]);
                        MatMulVecM8(MB_ij[2][j].Mat, temp_u8, &temp0[2]);
                        MatMulVecM8(MB_ij[3][j].Mat, temp_u8, &temp0[3]);

                        MatMulVecM8(Eij[r][0].Mat, temp0[0], &temp1[0]);
                        MatMulVecM8(Eij[r][1].Mat, temp0[1], &temp1[1]);
                        MatMulVecM8(Eij[r][2].Mat, temp0[2], &temp1[2]);
                        MatMulVecM8(Eij[r][3].Mat, temp0[3], &temp1[3]);
                        temp1[0].V ^= Eij[r][0].Vec.V;
                        temp1[1].V ^= Eij[r][1].Vec.V;
                        temp1[2].V ^= Eij[r][2].Vec.V;
                        temp1[3].V ^= Eij[r][3].Vec.V;
                        uint32_t m1 = temp1[0].V << 24 |
                                      temp1[1].V << 16 |
                                      temp1[2].V << 8 |
                                      temp1[3].V;
                        after_part2_change_linear[r][i][j][x][y]= 
                            (Out_before_xor[r][i][0][m1 >> 28 & 0xf] << 28) |
                            (Out_before_xor[r][i][1][m1 >> 24 & 0xf] << 24) |
                            (Out_before_xor[r][i][2][m1 >> 20 & 0xf] << 20) |
                            (Out_before_xor[r][i][3][m1 >> 16 & 0xf] << 16) |
                            (Out_before_xor[r][i][4][m1 >> 12 & 0xf] << 12) |
                            (Out_before_xor[r][i][5][m1 >> 8 & 0xf] << 8) |
                            (Out_before_xor[r][i][6][m1 >> 4 & 0xf] << 4) |
                            (Out_before_xor[r][i][7][m1 >> 0 & 0xf] << 0);
                    }
                }
            }
        }
    }

}