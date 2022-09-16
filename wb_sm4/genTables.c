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
            Table_addIn_part2[0][i][x] = x;

    for (int r = 1; r < 32; r++) {                                //进入s盒之前，去掉part1的非线性编码
        for(int i = 0; i < 4; i++){
            for (int x = 0; x < 256; x++) {
                uint8_t y = x;
                Table_addIn_part2[r][i][x] =
                    (In_part1[r][2 * i + 0][y >> 4 & 0xf] << 4) |
                    (In_part1[r][2 * i + 1][y >> 0 & 0xf] << 0);

            }
        }
    }


    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
                uint8_t temp_u8 = SBOX[Table_addIn_part2[i][j][x] ^ (ctx.sk[i] >> (24 - j * 8)) & 0xFF];   //异或rk，sbox
                uint32_t temp_u32 = temp_u8 << (24 - j * 8);                      //8bit扩展成32bit
                Table_SL[i][j][x] = MatMulNumM32(L_matrix, temp_u32);             //异或rk，sbox，L移位
            }
        }
    }

    for (int i = 0; i < 32; i++) {
        for (int j = 0; j < 4; j++) {
            for (int x = 0; x < 256; x++) {
               
                uint32_t temp_u32 = VecAddVecV32(Table_SL[i][j][x], Mask[i][j]);      //异或rk，sbox，L移位,+mask
                temp_u32 = MatMulNumM32(Ei[i].Mat, temp_u32);    //异或rk，sbox，L移位,+mask,+线性矩阵MB.Mat
                Table_part2[i][j][x] = VecAddVecV32(temp_u32, Ei[i].Vec.V);    //异或rk，sbox，L移位,+mask,+线性矩阵MB.Mat,+MB.Vec

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

}