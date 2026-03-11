# OpenMP parallelization added
#
# ビルド方法:
#   make          -> コンパイル
#   make run      -> dambreakサンプルを実行
#   make clean    -> オブジェクトファイルとバイナリを削除
#
# スレッド数の制御:
#   OMP_NUM_THREADS=4 ./mps_sim examples/dambreak/cal.txt
#
# コンパイラの自動検出:
#   - GCC (Linux 等) : -fopenmp のみで動作
#   - Apple Clang (macOS) : libomp が必要 (brew install libomp)
#
# Apple Clang 用の設定:
#   libomp が /opt/homebrew/opt/libomp にインストールされている前提。
#   インストール: brew install libomp

# --- コンパイラ・フラグの設定 ---
# macOS (Apple Clang) かどうかを検出して OpenMP フラグを切り替える
UNAME_S := $(shell uname -s)
CC       := gcc

ifeq ($(UNAME_S), Darwin)
    # macOS: Apple Clang は -fopenmp を直接サポートしないため libomp を使用
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(LIBOMP_PREFIX),)
        OMP_CFLAGS  := -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
        OMP_LDFLAGS := -L$(LIBOMP_PREFIX)/lib -lomp
    else
        $(warning libomp not found. Run: brew install libomp)
        OMP_CFLAGS  :=
        OMP_LDFLAGS :=
    endif
else
    # Linux / GCC
    OMP_CFLAGS  := -fopenmp
    OMP_LDFLAGS := -fopenmp
endif

CFLAGS  = -Wall -Wextra -O2 -I./include $(OMP_CFLAGS)
LDFLAGS = -lm $(OMP_LDFLAGS)

SRC_DIR = src
OBJ_DIR = obj
BIN     = mps_sim

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(SRCS:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

.PHONY: all clean run

all: $(BIN)

$(BIN): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

run: $(BIN)
	./$(BIN) examples/dambreak/cal.txt

clean:
	rm -rf $(OBJ_DIR) $(BIN)
