main.bin main.map main.rom: main.c util.c memory.c stack.asm integrals.c
	zcc \
	+embedded -clib=sdcc_iy \
	main.c terminal.c util.c util.asm memory.c integrals.c \
	-startup=1 \
	-pragma-define:CRT_ORG_CODE=0x1000 \
	-pragma-define:CRT_ORG_DATA=0x6100 \
	-pragma-define:REGISTER_SP=0x9FFF \
	-pragma-define:CRT_STACK_SIZE=256 \
	-pragma-define:CRT_INCLUDE_PREAMBLE=1 \
	-pragma-define:CLIB_FOPEN_MAX=0 \
	--max-allocs-per-node2000 \
	--opt-code-size \
	-SO3 -bn P2KHF.BIN \
	-create-app -lm \
	&& wc -c P2KHF.BIN \
	&& truncate -s 16k P2KHF.BIN