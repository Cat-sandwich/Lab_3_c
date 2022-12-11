#pragma once
#include <iostream>

class Exception
{
public:
	Exception();
	virtual void print() const = 0;
};
class Invalid_Index : public Exception
{
public:
	Invalid_Index();
	void print() const override;
};

class Empty : public Exception
{
public:
	Empty();
	void print() const override;
};

class Different_Dimensions : public Exception
{
public:
	Different_Dimensions();
	void print() const override;
};

class Divizion_By_Zero : public Exception
{
public:
	Divizion_By_Zero();
	void print() const override;
};

class Zero_Determinant : public Exception
{
public:
	Zero_Determinant();
	void print() const override;
};

class Dimensions_Incorrect : public Exception
{
public:
	Dimensions_Incorrect();
	void print() const override;
};


