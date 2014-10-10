/*
* Matt McCutchen's Big Integer Library
*/

/*
* This mechanism prevents files from being included twice.
* Each file gets its own `id' (here `NUMBERLIKEARRAY').
* When `#include'd, this file checks whether its `id' has
* already been flagged.  If not, it flags the `id' and
* loads the declarations.
*/
#ifndef NUMBERLIKEARRAY
#define NUMBERLIKEARRAY

// An essential memory-management constant.
// I wish this were built into C++ just as it is in Java.
#ifndef NULL
#define NULL 0
#endif

/*
* A NumberlikeArray<Blk> object holds a dynamically
* allocated array of Blk.  It provides certain basic
* memory management features needed by both BigUnsigned
* and BigUnsignedInABase, which are both derived from it.
*
* NumberlikeArray provides no information hiding, so make
* sure you know what you are doing if you use it directly.
* Classes derived from it will probably wish to pass on
* some members of NumberlikeArray to their clients while
* keeping some safe for themselves.  These classes should
* use protected inheritance and manually make some members
* public with declarations like this:
*
* public:
*     NumberlikeArray< whatever >::getLength;
*/

template <class Blk>
class NumberlikeArray {
	public:

	typedef unsigned int Index; // Type for the index of a block in the array
	static const unsigned int N; // The number of bits in a block, defined below.

	// FIELDS
	Index cap; // The current allocated capacity of this NumberlikeArray (in blocks)
	Index len; // The actual length of the value stored in this NumberlikeArray (in blocks)
	Blk *blk; // Dynamically allocated array of the blocks

	/*
	* Change made on 2005.01.06:
	*
	* If a zero-length NumberlikeArray is desired, no array is actually allocated.
	* Instead, `blk' is set to `NULL', and `cap' and `len' are zero as usual.
	*
	* `blk' is never dereferenced if the array has zero length.  Furthermore,
	* `delete NULL;' does nothing and causes no error. Therefore, we can use
	* `NULL' as if it were a zero-length array from `new'.
	*
	* This is a great convenience because the only code that need be changed
	* is the array allocation code.  All other code will still work fine.
	*/

	// MANAGEMENT
	NumberlikeArray(Index c) : cap(c), len(0) { // Creates a NumberlikeArray with a capacity
		blk = (cap > 0) ? (new Blk[cap]) : NULL;
	}
	void allocate(Index c); // Ensures the array has at least the indicated capacity, maybe discarding contents
	void allocateAndCopy(Index c); // Ensures the array has at least the indicated capacity, preserving its contents

	/*
	* Default constructor.
	*
	* If a class derived from NumberlikeArray knows at initializer time what size array
	* it wants, it can call the first constructor listed above in an initializer.
	*
	* Otherwise, this default constructor will be implicitly invoked, pointing `blk' to
	* `NULL', a fake zero-length block array.  The derived class can allocate the desired
	* array itself and overwrite `blk'; it need not `delete [] blk' first.
	*
	* This change fixes a memory leak reported by Milan Tomic on 2005.01.06.
	* Integer-type-to-BigUnsigned (and BigInteger) conversion constructors have always
	* allocated their own array of length 0 or 1 after seeing whether the input is zero.
	* But when the NumberlikeArray transition occurred, these constructors contained an
	* implicit initializer call to the old NumberlikeArray default constructor, which
	* created a real `new'-allocated zero-length array.  This array would then be lost,
	* causing a small but annoying memory leak.
	*/
	NumberlikeArray() : cap(0), len(0) {
		blk = NULL;
	}
	NumberlikeArray(const NumberlikeArray<Blk> &x); // Copy constructor
	void operator=(const NumberlikeArray<Blk> &x); // Assignment operator
	NumberlikeArray(const Blk *b, Index l); // Constructor from an array of blocks
	~NumberlikeArray() { // Destructor
		delete [] blk; // Does nothing and causes no error if `blk' is null.
	}

	// PICKING APART
	// These accessors can be used to get the pieces of the value
	Index getCapacity() const { return cap; }
	Index getLength() const { return len; }
	Blk getBlock(Index i) const { return blk[i]; };
	bool isEmpty() const { return len == 0; }

	// Equality comparison: checks if arrays have same length and matching values
	// Derived classes may wish to override these if differing arrays can
	// sometimes be considered equivalent.
	bool operator ==(const NumberlikeArray<Blk> &x) const;
	bool operator !=(const NumberlikeArray<Blk> &x) const { return !operator ==(x);	}

};

/*
* =================================
* BELOW THIS POINT are template definitions; above are declarations.
*
* Definitions would ordinarily belong in a file NumberlikeArray.cc so that they would
* be compiled once into NumberlikeArray.o and then linked.
*
* However, because of the way templates are usually implemented,
* template ``definitions'' are treated as declarations by the compiler.
* When someone uses an instance of the template, definitions are generated,
* and the linker is smart enough to toss duplicate definitions for the same
* instance generated by different files.
*
* Thus, the template ``definitions'' for NumberlikeArray must appear in this header file
* so other files including NumberlikeArray will be able to generate real definitions.
*/

template <class Blk>
const unsigned int NumberlikeArray<Blk>::N = 8 * sizeof(Blk);

// MANAGEMENT

// This routine is called to ensure the array is at least a
// certain size before another value is written into it.
template <class Blk>
void NumberlikeArray<Blk>::allocate(Index c) {
	// If the requested capacity is more than the current capacity...
	if (c > cap) {
		// Delete the old number array
		delete [] blk;
		// Allocate the new array
		cap = c;
		blk = new Blk[cap];
	}
}

// This routine is called to ensure the array is at least a
// certain size without losing its contents.
template <class Blk>
void NumberlikeArray<Blk>::allocateAndCopy(Index c) {
	// If the requested capacity is more than the current capacity...
	if (c > cap) {
		Blk *oldBlk = blk;
		// Allocate the new number array
		cap = c;
		blk = new Blk[cap];
		// Copy number blocks
		Index i;
		for (i = 0; i < len; i++)
			blk[i] = oldBlk[i];
		// Delete the old array
		delete [] oldBlk;
	}
}

// Copy constructor
template <class Blk>
NumberlikeArray<Blk>::NumberlikeArray(const NumberlikeArray<Blk> &x) : len(x.len) {
	// Create array
	cap = len;
	blk = new Blk[cap];
	// Copy blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = x.blk[i];
}

// Assignment operator
template <class Blk>
void NumberlikeArray<Blk>::operator=(const NumberlikeArray<Blk> &x) {
	// Calls like a = a have no effect
	if (this == &x)
		return;
	// Copy length
	len = x.len;
	// Expand array if necessary
	allocate(len);
	// Copy number blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = x.blk[i];
}

// Constructor from an array of blocks
template <class Blk>
NumberlikeArray<Blk>::NumberlikeArray(const Blk *b, Index l) : cap(l), len(l) {
	// Create array
	blk = new Blk[cap];
	// Copy blocks
	Index i;
	for (i = 0; i < len; i++)
		blk[i] = b[i];
}


// EQUALITY TEST
// This uses == to compare Blks for equality.
// Therefore, Blks must have an == operator with the desired semantics.
template <class Blk>
bool NumberlikeArray<Blk>::operator ==(const NumberlikeArray<Blk> &x) const {
	// Different lengths imply different objects.
	if (len != x.len)
		return false;
	else {
		// Compare matching blocks one by one.
		Index i;
		for (i = 0; i < len; i++)
			if (blk[i] != x.blk[i])
				return false;
		// If no blocks differed, the objects are equal.
		return true;
	}
}

#endif
/*
* Matt McCutchen's Big Integer Library
*/

#ifndef BIGUNSIGNED
#define BIGUNSIGNED

//#include "NumberlikeArray.hh"

/*
* A BigUnsigned object represents a nonnegative integer of size
* limited only by available memory.  A BigUnsigned can be
* created from and converted back to most integral types,
* and many math operations are defined on BigUnsigneds.
*
* The number is stored as a series of blocks in a
* dynamically allocated array.  It is as if the number
* were written digit by digit in base 2 ^ N, **where N is the
* number of bits in an unsigned long.**
*
* The memory-management details that used to be in here have
* been moved into NumberlikeArray, which BigUnsigned now derives from.
* `(NlA)' means that member(s) are declared identically in NumberlikeArray.
* Such members are either redeclared here to make them public or are
* here, commented out, for reference.
*/

class BigUnsigned : protected NumberlikeArray<unsigned long> {

	// TYPES & CONSTANTS
	public:
	enum CmpRes { less = -1, equal = 0, greater = 1 }; // Enumeration for the result of a comparison
	typedef unsigned long Blk; // The number block type that BigUnsigneds are built from
	typedef NumberlikeArray<Blk>::Index Index; // (NlA) Type for the index of a block in the array
	NumberlikeArray<Blk>::N; // Number of bits in a Blk

	/*
	// FIELDS
	protected:
	Index cap; // (NlA) The current allocated capacity of this BigUnsigned (in blocks)
	Index len; // (NlA) The actual length of the number stored in this BigUnsigned (in blocks)
	Blk *blk; // (NlA) Dynamically allocated array of the number blocks
	*/

	// MANAGEMENT
	protected:
	// These members generally defer to those in NumberlikeArray, possibly with slight changes.
	// It might be nice if one could request that constructors be inherited in C++.

	BigUnsigned(int, Index c) : NumberlikeArray<Blk>(0, c) {} // Creates a BigUnsigned with a capacity

	void zapLeadingZeros() { // Decreases len to eliminate leading zeros
		while (len > 0 && blk[len - 1] == 0)
			len--;
	}

	//void allocate(Index c); // (NlA) Ensures the number array has at least the indicated capacity, maybe discarding contents
	//void allocateAndCopy(Index c); // (NlA) Ensures the number array has at least the indicated capacity, preserving its contents

	public:
	BigUnsigned() : NumberlikeArray<Blk>() {} // Default constructor (value is 0)
	BigUnsigned(const BigUnsigned &x) : NumberlikeArray<Blk>(x) {} // Copy constructor

	void operator=(const BigUnsigned &x) { // Assignment operator
		NumberlikeArray<Blk>::operator =(x);
	}

	BigUnsigned(const Blk *b, Index l) : NumberlikeArray<Blk>(b, l) { // Constructor from an array of blocks
		zapLeadingZeros();
	}

	// Constructors from integral types
	BigUnsigned(unsigned long  x);
	BigUnsigned(         long  x);
	BigUnsigned(unsigned int   x);
	BigUnsigned(         int   x);
	BigUnsigned(unsigned short x);
	BigUnsigned(         short x);
	~BigUnsigned() {} // Destructor

	// CONVERTERS to integral types
	public:
	operator unsigned long () const;
	operator          long () const;
	operator unsigned int  () const;
	operator          int  () const;
	operator unsigned short() const;
	operator          short() const;

	// PICKING APART
	// These accessors can be used to get the pieces of the number
	public:
	NumberlikeArray<Blk>::getCapacity;
	NumberlikeArray<Blk>::getLength;
	// Note that getBlock returns 0 if the block index is beyond the length of the number.
	// A routine that uses this accessor can safely assume a BigUnsigned has 0s infinitely to the left.
	Blk getBlock(Index i) const { return i >= len ? 0 : blk[i]; }
	// Note how we replace one level of abstraction with another.  Isn't that neat?
	bool isZero() const { return NumberlikeArray<Blk>::isEmpty(); } // Often convenient for loops

	// COMPARISONS
	public:
	// Compares this to x like Perl's <=>
	CmpRes compareTo(const BigUnsigned &x) const;
	// Normal comparison operators
	// Bug fixed 2006.04.24: Only we, not the user, can pass a BigUnsigned off as a
	// NumberlikeArray, so we have to wrap == and !=.
	bool operator ==(const BigUnsigned &x) const {
		return NumberlikeArray<Blk>::operator ==(x);
	}
	bool operator !=(const BigUnsigned &x) const {
		return NumberlikeArray<Blk>::operator !=(x);
	}
	bool operator < (const BigUnsigned &x) const { return compareTo(x) == less   ; }
	bool operator <=(const BigUnsigned &x) const { return compareTo(x) != greater; }
	bool operator >=(const BigUnsigned &x) const { return compareTo(x) != less   ; }
	bool operator > (const BigUnsigned &x) const { return compareTo(x) == greater; }

	/*
	* BigUnsigned and BigInteger both provide three kinds of operators.
	* Here ``big-integer'' refers to BigInteger or BigUnsigned.
	*
	* (1) Overloaded ``return-by-value'' operators:
	*     +, -, *, /, %, unary -.
	* Big-integer code using these operators looks identical to
	* code using the primitive integer types.  These operators take
	* one or two big-integer inputs and return a big-integer result,
	* which can then be assigned to a BigInteger variable or used
	* in an expression.  Example:
	*     BigInteger a(1), b = 1;
	*     BigInteger c = a + b;
	*
	* (2) Overloaded assignment operators:
	*     +=, -=, *=, /=, %=, &=, |=, ^=, ++, --, flipSign.
	* Again, these are used on big integers just like on ints.
	* They take one writable big integer that both provides an
	* operand and receives a result.  The first eight also take
	* a second read-only operand.  Example:
	*     BigInteger a(1), b(1);
	*     a += b;
	*
	* (3) ``Put-here'' operations: `add', `subtract', etc.
	* Using a return-by-value or assignment operator generally involves
	* copy constructions and/or assignments.  The ``put-here'' operations
	* require none, but they are more of a hassle to use.  Most take two
	* read-only operands and save the result in the calling object `*this',
	* whose previous value is ignored.  `divideWithRemainder' is an exception.
	* <<< NOTE >>>: Put-here operations do not return a value: they don't need to!!
	* Examples:
	*     BigInteger a(43), b(7), c, d;
	*     c = a + b;   // Now c == 50.
	*     c.add(a, b); // Same effect but without the two bulk-copies.
	*     c.divideWithRemainder(b, d); // 50 / 7; now d == 7 (quotient) and c == 1 (remainder).
	*     a.add(a, b); // ``Aliased'' calls now do the right thing using a
	*              // temporary copy, but see note on divideWithRemainder.
	*/

	// PUT-HERE OPERATIONS
	public:
	/* These 3: Two read-only operands as arguments.  Result left in *this. */
	void add(const BigUnsigned &a, const BigUnsigned &b); // Addition
	void subtract(const BigUnsigned &a, const BigUnsigned &b); // Subtraction
	void multiply(const BigUnsigned &a, const BigUnsigned &b); // Multiplication
	/* Divisive stuff
	* `a.divideWithRemainder(b, q)' is like `q = a / b, a %= b'.
	* Semantics similar to Donald E. Knuth's are used for / and %,
	* and these differ from the semantics of primitive-type
	* / and % under division by zero.
	* Look in `BigUnsigned.cc' for details.
	* `a.divideWithRemainder(b, a)' causes an exception: it doesn't make
	* sense to write quotient and remainder into the same variable.
	*/
	void divideWithRemainder(const BigUnsigned &b, BigUnsigned &q);
	void divide(const BigUnsigned &a, const BigUnsigned &b) {
		BigUnsigned a2(a);
		a2.divideWithRemainder(b, *this);
		// quotient now in *this
		// don't care about remainder left in a2
	}
	void modulo(const BigUnsigned &a, const BigUnsigned &b) {
		*this = a;
		BigUnsigned q;
		divideWithRemainder(b, q);
		// remainder now in *this
		// don't care about quotient left in q
	}
	// Bitwise operations.  Result left in *this.
	// These are not provided for BigIntegers; I think that using them on BigIntegers
	// will discard the sign first.
	void bitAnd(const BigUnsigned &a, const BigUnsigned &b); // Bitwise AND
	void bitOr(const BigUnsigned &a, const BigUnsigned &b); // Bitwise OR
	void bitXor(const BigUnsigned &a, const BigUnsigned &b); // Bitwise XOR
	void bitShiftLeft(const BigUnsigned &a, unsigned int b); // Bitwise left shift
	void bitShiftRight(const BigUnsigned &a, unsigned int b); // Bitwise right shift

	// NORMAL OPERATORS
	// These perform the operation on this (to the left of the operator)
	// and x (to the right of the operator) and return a new BigUnsigned with the result.
	public:
	BigUnsigned operator +(const BigUnsigned &x) const; // Addition
	BigUnsigned operator -(const BigUnsigned &x) const; // Subtraction
	BigUnsigned operator *(const BigUnsigned &x) const; // Multiplication
	BigUnsigned operator /(const BigUnsigned &x) const; // Division
	BigUnsigned operator %(const BigUnsigned &x) const; // Modular reduction
	BigUnsigned operator &(const BigUnsigned &x) const; // Bitwise AND
	BigUnsigned operator |(const BigUnsigned &x) const; // Bitwise OR
	BigUnsigned operator ^(const BigUnsigned &x) const; // Bitwise XOR
	BigUnsigned operator <<(unsigned int b) const; // Bitwise left shift
	BigUnsigned operator >>(unsigned int b) const; // Bitwise right shift
	// Additional operators in an attempt to avoid overloading tangles.
	BigUnsigned operator <<(int b) const;
	BigUnsigned operator >>(int b) const;

	// ASSIGNMENT OPERATORS
	// These perform the operation on this and x, storing the result into this.
	public:
	void operator +=(const BigUnsigned &x); // Addition
	void operator -=(const BigUnsigned &x); // Subtraction
	void operator *=(const BigUnsigned &x); // Multiplication
	void operator /=(const BigUnsigned &x); // Division
	void operator %=(const BigUnsigned &x); // Modular reduction
	void operator &=(const BigUnsigned &x); // Bitwise AND
	void operator |=(const BigUnsigned &x); // Bitwise OR
	void operator ^=(const BigUnsigned &x); // Bitwise XOR
	void operator <<=(unsigned int b); // Bitwise left shift
	void operator >>=(unsigned int b); // Bitwise right shift
	// Additional operators in an attempt to avoid overloading tangles.
	void operator <<=(int b);
	void operator >>=(int b);

	// INCREMENT/DECREMENT OPERATORS
	// These increase or decrease the number by 1.  To discourage side effects,
	// these do not return *this, so prefix and postfix behave the same.
	public:
	void operator ++(   ); // Prefix  increment
	void operator ++(int); // Postfix decrement
	void operator --(   ); // Prefix  increment
	void operator --(int); // Postfix decrement

	// Helper function that needs access to BigUnsigned internals
	friend Blk getShiftedBlock(const BigUnsigned &num, Index x, unsigned int y);
};

// NORMAL OPERATORS
/* These create an object to hold the result and invoke
* the appropriate put-here operation on it, passing
* this and x.  The new object is then returned. */
inline BigUnsigned BigUnsigned::operator +(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.add(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator -(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.subtract(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator *(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.multiply(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator /(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.divide(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator %(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.modulo(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator &(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.bitAnd(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator |(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.bitOr(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator ^(const BigUnsigned &x) const {
	BigUnsigned ans;
	ans.bitXor(*this, x);
	return ans;
}
inline BigUnsigned BigUnsigned::operator <<(unsigned int b) const {
	BigUnsigned ans;
	ans.bitShiftLeft(*this, b);
	return ans;
}
inline BigUnsigned BigUnsigned::operator >>(unsigned int b) const {
	BigUnsigned ans;
	ans.bitShiftRight(*this, b);
	return ans;
}
inline BigUnsigned BigUnsigned::operator <<(int b) const {
	if (b < 0)
		throw "BigUnsigned::operator <<(int): Negative shift amounts are not supported";
	return *this << (unsigned int)(b);
}
inline BigUnsigned BigUnsigned::operator >>(int b) const {
	if (b < 0)
		throw "BigUnsigned::operator >>(int): Negative shift amounts are not supported";
	return *this >> (unsigned int)(b);
}

/*
 * ASSIGNMENT OPERATORS
 *
 * Now the responsibility for making a temporary copy if necessary
 * belongs to the put-here operations.  I made this change on 2007.02.13 after
 * Boris Dessy pointed out that the old implementation handled calls like
 * "a *= a" badly: it translated them to essentially "a.multiply(aCopy, a)",
 * which threw an exception.
 */
inline void BigUnsigned::operator +=(const BigUnsigned &x) {
	add(*this, x);
}
inline void BigUnsigned::operator -=(const BigUnsigned &x) {
	subtract(*this, x);
}
inline void BigUnsigned::operator *=(const BigUnsigned &x) {
	multiply(*this, x);
}
inline void BigUnsigned::operator /=(const BigUnsigned &x) {
	// Updated for divideWithRemainder
	BigUnsigned thisCopy(*this);
	thisCopy.divideWithRemainder(x, *this);
	// quotient left in *this
	// don't care about remainder left in thisCopy
}
inline void BigUnsigned::operator %=(const BigUnsigned &x) {
	// Shortcut (woohoo!)
	BigUnsigned q;
	divideWithRemainder(x, q);
	// remainder left in *this
	// don't care about quotient left in q
}
inline void BigUnsigned::operator &=(const BigUnsigned &x) {
	bitAnd(*this, x);
}
inline void BigUnsigned::operator |=(const BigUnsigned &x) {
	bitOr(*this, x);
}
inline void BigUnsigned::operator ^=(const BigUnsigned &x) {
	bitXor(*this, x);
}
inline void BigUnsigned::operator <<=(unsigned int b) {
	bitShiftLeft(*this, b);
}
inline void BigUnsigned::operator >>=(unsigned int b) {
	bitShiftRight(*this, b);
}
inline void BigUnsigned::operator <<=(int b) {
	if (b < 0)
		throw "BigUnsigned::operator <<=(int): Negative shift amounts are not supported";
	*this <<= (unsigned int)(b);
}
inline void BigUnsigned::operator >>=(int b) {
	if (b < 0)
		throw "BigUnsigned::operator >>=(int): Negative shift amounts are not supported";
	*this >>= (unsigned int)(b);
}

#endif
/*
* Matt McCutchen's Big Integer Library
*/

#ifndef BIGUNSIGNEDINABASE
#define BIGUNSIGNEDINABASE

//#include "NumberlikeArray.hh"
//#include "BigUnsigned.hh"
#include <string>

/*
* A BigUnsignedInABase object represents a nonnegative
* integer of size limited only by available memory,
* represented in a user-specified base that can fit in
* an `unsigned short' (most can, and this saves memory).
*
* BigUnsignedInABase is intended as an intermediary class
* with little functionality of its own.  BigUnsignedInABase
* objects can be constructed from, and converted to,
* BigUnsigneds (requiring multiplication, mods, etc.) and
* `std::string's (by switching digit values for appropriate
* characters).
*
* BigUnsignedInABase is similar to BigUnsigned.  Note the following:
*
* (1) They represent the number in exactly the same way, except
* that BigUnsignedInABase uses ``digits'' (or Digit) where BigUnsigned uses
* ``blocks'' (or Blk).
*
* (2) Both use the management features of NumberlikeArray.  (In fact,
* my desire to add a BigUnsignedInABase class without duplicating a
* lot of code led me to introduce NumberlikeArray.)
*
* (3) The only arithmetic operation supported by BigUnsignedInABase
* is an equality test.  Use BigUnsigned for arithmetic.
*/

class BigUnsignedInABase : protected NumberlikeArray<unsigned short> {

	// TYPES
	public:
	typedef unsigned short Digit; // The digit type that BigUnsignedInABases are built from
	typedef Digit Base;

	// FIELDS
	protected:
	Base base; // The base of this BigUnsignedInABase

	// MANAGEMENT
	protected:
	// These members generally defer to those in NumberlikeArray, possibly with slight changes.
	// It might be nice if one could request that constructors be inherited in C++.

	BigUnsignedInABase(int, Index c) : NumberlikeArray<Digit>(0, c) {} // Creates a BigUnsignedInABase with a capacity

	void zapLeadingZeros() { // Decreases len to eliminate leading zeros
		while (len > 0 && blk[len - 1] == 0)
			len--;
	}

	//void allocate(Index c); // (NlA) Ensures the number array has at least the indicated capacity, maybe discarding contents
	//void allocateAndCopy(Index c); // (NlA) Ensures the number array has at least the indicated capacity, preserving its contents

	public:
	BigUnsignedInABase() : NumberlikeArray<Digit>(), base(2) {} // Default constructor (value is 0 in base 2)
	BigUnsignedInABase(const BigUnsignedInABase &x) : NumberlikeArray<Digit>(x), base(x.base) {} // Copy constructor

	void operator =(const BigUnsignedInABase &x) { // Assignment operator
		NumberlikeArray<Digit>::operator =(x);
		base = x.base;
	}

	BigUnsignedInABase(const Digit *d, Index l) : NumberlikeArray<Digit>(d, l) { // Constructor from an array of digits
		zapLeadingZeros();
	}

	// LINKS TO BIGUNSIGNED
	BigUnsignedInABase(const BigUnsigned &x, Base base);
	operator BigUnsigned() const;

	/* LINKS TO STRINGS
	*
	* These use the symbols ``0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'' to represent
	* digits of 0 through 35.  When parsing strings, lowercase is also accepted.
	*
	* All string representations are big-endian (big-place-value digits first).
	* (Computer scientists have adopted zero-based counting; why can't they
	* tolerate little-endian numbers?  It makes a lot of sense!)
	*
	* No string representation has a ``base indicator'' like ``0x''.
	*
	* An exception is made for zero: it is converted to ``0'' and not the empty string.
	*
	* If you want different conventions, write your
	* own routines to go between BigUnsignedInABase and strings.  It's not hard.
	*/
	operator std::string() const;
	BigUnsignedInABase(const std::string &s, Base base);

	// PICKING APART
	// These accessors can be used to get the pieces of the number
	public:
	Base getBase() const { return base; }
	NumberlikeArray<Digit>::getCapacity; // (NlA)
	NumberlikeArray<Digit>::getLength; // (NlA)
	// Note that getDigit returns 0 if the digit index is beyond the length of the number.
	// A routine that uses this accessor can safely assume a BigUnsigned has 0s infinitely to the left.
	Digit getDigit(Index i) const { return i >= len ? 0 : blk[i]; }
	// Note how we replace one level of abstraction with another.
	bool isZero() const { return NumberlikeArray<Digit>::isEmpty(); } // Often convenient for loops

	// EQUALITY TEST
	public:
	// Equality test
	bool operator ==(const BigUnsignedInABase &x) const {
		return base == x.base && NumberlikeArray<Digit>::operator ==(x);
	}
	bool operator !=(const BigUnsignedInABase &x) const { return !operator ==(x); }

};

#endif

