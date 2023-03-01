# 符号ベース暗号に対する安全性解析

本プロジェクトでは，令和５年千葉大学大学院　融合理工学府における修士論文「符号ベース暗号に対する安全性解析」7.3節，7.4節で用いたコードを提供する．

## ディレクトリ構成

```
  .
  ├── Hamming_CBC
  │   ├── Hamming_CBC_params.c
  │   ├── Hamming_CBC_params.h
  │   ├── Hamming_CBC_params.o
  │   ├── MMT_BJMM
  │   │   ├── Dicke_quantum_MMT_BJMM
  │   │   ├── Dicke_quantum_MMT_BJMM.c
  │   │   ├── Dicke_quantum_MMT_BJMM.o
  │   │   ├── MMT_BJMM.h
  │   │   ├── MMT_BJMM_params.c
  │   │   ├── MMT_BJMM_params.h
  │   │   ├── MMT_BJMM_params.o
  │   │   ├── classical_MMT_BJMM
  │   │   ├── classical_MMT_BJMM.c
  │   │   ├── classical_MMT_BJMM.o
  │   │   ├── improved_quantum_MMT_BJMM
  │   │   ├── improved_quantum_MMT_BJMM.c
  │   │   ├── improved_quantum_MMT_BJMM.o
  │   │   ├── quantum_MMT_BJMM
  │   │   ├── quantum_MMT_BJMM.c
  │   │   ├── quantum_MMT_BJMM.o
  │   │   ├── quantum_MMT_BJMM_classical_PRAM
  │   │   ├── quantum_MMT_BJMM_classical_PRAM.c
  │   │   ├── quantum_MMT_BJMM_classical_PRAM.o
  │   │   ├── quantum_MMT_BJMM_classical_PRAM_Tcost
  │   │   ├── quantum_MMT_BJMM_classical_PRAM_Tcost.c
  │   │   └── quantum_MMT_BJMM_classical_PRAM_Tcost.o
  │   └── Prange
  │       ├── Prange.h
  │       ├── Prange_params.c
  │       ├── Prange_params.h
  │       ├── Prange_params.o
  │       ├── classical_Prange
  │       ├── classical_Prange.c
  │       ├── classical_Prange.o
  │       ├── improved_quantum_Prange
  │       ├── improved_quantum_Prange.c
  │       ├── improved_quantum_Prange.o
  │       ├── quantum_Prange
  │       ├── quantum_Prange.c
  │       ├── quantum_Prange.o
  │       ├── quantum_Prange_classical_PRAM
  │       ├── quantum_Prange_classical_PRAM.c
  │       ├── quantum_Prange_classical_PRAM.o
  │       ├── quantum_Prange_classical_PRAM_Tcost
  │       ├── quantum_Prange_classical_PRAM_Tcost.c
  │       └── quantum_Prange_classical_PRAM_Tcost.o
  ├── Rank_CBC
  │   ├── GRS
  │   │   ├── GRS.h
  │   │   ├── GRS_params.c
  │   │   ├── GRS_params.h
  │   │   ├── GRS_params.o
  │   │   ├── classical_GRS
  │   │   ├── classical_GRS.c
  │   │   ├── classical_GRS.o
  │   │   ├── improved_quantum_GRS
  │   │   ├── improved_quantum_GRS.c
  │   │   ├── improved_quantum_GRS.o
  │   │   ├── quantum_GRS
  │   │   ├── quantum_GRS.c
  │   │   ├── quantum_GRS.o
  │   │   ├── quantum_GRS_classical_PRAM
  │   │   ├── quantum_GRS_classical_PRAM.c
  │   │   └── quantum_GRS_classical_PRAM.o
  │   ├── Rank_CBC_params.c
  │   ├── Rank_CBC_params.h
  │   └── Rank_CBC_params.o
  ├── basic_operation.c
  ├── basic_operation.h
  └── basic_operation.o
```

------------------------------

## 7.3節のコードの概説
7.3節では，10個の攻撃アルゴリズムによる計算コストの導出を議論した．それぞれの項で説明されるアルゴリズムがどのコードに対応するかを説明する．

| 項番号  | 項タイトル | 対応するコード |
| ------ | ------- | ----------- |
| 7.3.1  | 古典回路上での量子 Prange アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/classical_Prange* |
| 7.3.2  | 量子回路上での量子 Prange アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_Prange* |
| 7.3.3  | 量子回路上での改善版量子 Prange アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/improved_quantum_Prange* |
| 7.3.4  | 古典回路上での量子 MMT/BJMM アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/classical_MMT_BJMM* |
| 7.3.5  | Dicke 状態を用いない場合の量子回路上での量子 MMT/BJMM アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_MMT_BJMM* |
| 7.3.6  | Dicke 状態を用いる場合の量子回路上での量子 MMT/BJMM アルゴリズムの計算 | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/Dicke_quantum_MMT_BJMM* |
| 7.3.7  | 量子回路上での改善版量子 MMT/BJMM アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/improved_quantum_MMT_BJMM* |
| 7.3.8  | 古典回路上での量子 GRS アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/classical_GRS* |
| 7.3.9  | 量子回路上での量子 GRS アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_GRS* |
| 7.3.10 | 量子回路上での改善版量子 GRS アルゴリズムの計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/improved_quantum_GRS* |

-----------------------------------------

## 7.4節のコードの概説
7.4節では，いくつかの古典PRAM演算による計算コストを解説した．それぞれの計算コストがどのコードに対応するかを下記の表に示す．

| 計算コスト  | 対応するコード |
| ------ | ----------- |
| 量子回路上の量子 Prange アルゴリズムの古典 PRAM 演算による計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_Prange_classical_PRAM* |
| 量子回路上の量子 MMT/BJMM アルゴリズムの古典 PRAM 演算による計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_MMT_BJMM_classical_PRAM* |
| 量子回路上の量子 GRS アルゴリズムの古典 PRAM 演算による計算コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_GRS_classical_PRAM* |
| 量子回路上の量子 Prange アルゴリズムの古典 PRAM 演算による T ゲートベースの DW コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_Prange_classical_PRAM_Tcost* |
| 量子回路上の量子 MMT/BJMM アルゴリズムの古典 PRAM 演算による T ゲートベースの DW コスト | ./basic_operation.* <br> ./Hamming_CBC/Hamming_CBC_params.* <br> ./Hamming_CBC/Prange/quantum_MMT_BJMM_classical_PRAM_Tcost* |

