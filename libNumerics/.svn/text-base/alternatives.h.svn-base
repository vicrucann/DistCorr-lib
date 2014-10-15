#ifndef ALTERNATIVES_H
#define ALTERNATIVES_H

#ifdef _WIN32
#if (!_WIN64)
inline double ddot_asm(double* v1_ptr, double* v2_ptr, int cnt) {
	double res = 0;
	__asm
	{
		push esi 
		push edi

		mov ecx, cnt
		mov esi, v1_ptr
		mov edi, v2_ptr

		_loop:
			fld qword ptr [esi]
			fld qword ptr [edi]

			fmulp st(1), st
			fld qword ptr res
			faddp st(1), st
			fstp qword ptr res
		
			add esi, 8
			add edi, 8
			dec ecx
		jnz _loop
		
		pop edi
		pop esi
	}
	return res;
}
#endif
#endif

#endif
