# 符号ベース暗号に対する安全性解析

本プロジェクトでは，令和５年千葉大学大学院　融合理工学府における修士論文「符号ベース暗号に対する安全性解析」7.3節，7.4節で用いたコードを提供する．

## ディレクトリ構成

<pre>
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
 <pre>

------------------------------

## それぞれのコードの概説
7.3節では，10個の攻撃アルゴリズムによる計算コストの導出を議論した．それぞれの項で説明されるアルゴリズムが上記のコードのいずれに対応するかを説明する．

| 項番号  | 項タイトル | 対応するコード |
| ------ | ------- | ----------- |
| 7.3.1  | 古典回路上での量子 Prange アルゴリズムの計算コスト | ./basic_operation.*, ./Hamming_CBC/Hamming_CBC_params.*, ./Hamming_CBC/Prange/classical_Prange* |
