
#ifndef BITBOOLS_H
#define BITBOOLS_H

namespace xyy
{
	static unsigned char TrueMasks[8] = { 1, 2, 4, 8, 16, 32, 64, 128 };
	static unsigned char FalseMasks[8] = { 254, 253, 251, 247, 239, 223, 191, 127 };

	class BitBools
	{
	protected:
		unsigned char *_bits;
		int            _m;

	public:
		BitBools() : _bits(NULL), _m(0) {}
		BitBools(int n, bool initvalue) : _bits(NULL), _m(0)
		{
			set_number(n, initvalue);
		}
		BitBools(const BitBools &rhs)
		{
			_m = rhs._m;
			_bits = (unsigned char*)malloc(_m * sizeof(unsigned char));

			memcpy(_bits, rhs._bits, sizeof(unsigned char) * _m);
		}
		BitBools& operator= (const BitBools &rhs)
		{
			clear();

			_m = rhs._m;
			_bits = (unsigned char*)malloc(sizeof(unsigned char) * _m);

			memcpy(_bits, rhs._bits, sizeof(unsigned char) * _m);

			return *this;
		}

		~BitBools()
		{ 
			clear();
		}

		void clear()
		{
			if (_bits)
				free(_bits);

			_bits = NULL;
			_m = 0;
		}

		void set_number(int n, bool initvalue)
		{
			clear();

			_m = n / 8 + 1;
			_bits = (unsigned char*)calloc(_m, sizeof(unsigned char));

			if (initvalue)
				memset(_bits, 1, sizeof(unsigned char) * _m);
		}

		void reset(bool initvalue)
		{
			if (initvalue)
				memset(_bits, 1, sizeof(unsigned char) * _m);
			else
				memset(_bits, 0, sizeof(unsigned char) * _m);
		}
		
		bool is_true(int k)
		{
			return (_bits[k / 8] & TrueMasks[k % 8]) > 0;
		}

		void set_true(int k)
		{
			_bits[k / 8] |= TrueMasks[k % 8];
		}

		void set_false(int k)
		{
			_bits[k / 8] &= FalseMasks[k % 8];
		}
	};
}

#endif
