// //////////////////////////////////////////////////////////
// MemoryMapped.cpp
// Copyright (c) 2013 Stephan Brumme. All rights reserved.
// see http://create.stephan-brumme.com/disclaimer.html
// Modified by Daniel Siladji to use mman.h on Windows, as well

#include "MemoryMapped.h"

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>

#ifdef __unix
#define O_BINARY 0
#include <unistd.h>
#include <sys/mman.h>
#else
#include "mman.h"
#include <io.h>
#include <windows.h>
#endif

#ifdef _MSC_VER
//auto open = _open;
//auto creat = _creat;
auto lseek64 = _lseeki64;
auto fstat64 = _fstat64;
#define stat64 _stat64
#define O_CREAT _O_CREAT
#define O_BINARY _O_BINARY
#define O_RDWR _O_RDWR
#define O_TRUNC _O_TRUNC
#define S_IWRITE _S_IWRITE
#define S_IREAD _S_IREAD
#endif

/// do nothing, must use open()
MemoryMapped::MemoryMapped()
	: _filename(),
	_filesize(0),
	_hint(Normal),
	_mappedBytes(0),
	_file(0),
	_mappedView(NULL)
{
}


/// open file, mappedBytes = 0 maps the whole file
MemoryMapped::MemoryMapped(std::string_view filename, size_t mappedBytes, CacheHint hint)
	: _filename(filename),
	_filesize(0),
	_hint(hint),
	_mappedBytes(mappedBytes),
	_file(0),
	_mappedView(NULL)
{
	open(filename, mappedBytes, hint);
}


/// close file (see close() )
MemoryMapped::~MemoryMapped()
{
	close();
}


/// open file
bool MemoryMapped::open(std::string_view filename, size_t mappedBytes, CacheHint hint)
{
	auto mask = O_CREAT | O_BINARY | O_RDWR;
	if (mappedBytes != 0) mask |= O_TRUNC;

	_file = ::open(filename.data(), mask , S_IWRITE | S_IREAD);
	if (_file == -1)
	{
		_file = 0;
		return false;
	}

	// file size
	if (mappedBytes == 0) {
		struct stat64 statInfo;
		if (fstat64(_file, &statInfo) < 0)
			return false;

		_filesize = statInfo.st_size;
	}
	else {
		lseek64(_file, mappedBytes - 1, SEEK_SET);
		write(_file, "", 1);
		_filesize = mappedBytes;
	}
	
	// initial mapping
	remap(0, _filesize);

	if (!_mappedView)
		return false;

	// everything's fine
	return true;
}


/// close file
void MemoryMapped::close()
{
	munmap(_mappedView, _mappedBytes);
}


/// access position, no range checking (faster)
unsigned char& MemoryMapped::operator[](size_t offset)
{
	return ((unsigned char*)_mappedView)[offset];
}


/// access position, including range checking
unsigned char& MemoryMapped::at(size_t offset)
{
	// checks
	if (!_mappedView)
		throw std::invalid_argument("No view mapped");
	if (offset >= _filesize)
		throw std::out_of_range("View is not large enough");

	return operator[](offset);
}


/// raw access
unsigned char* MemoryMapped::getData() const
{
	return (unsigned char*)_mappedView;
}


/// true, if file successfully opened
bool MemoryMapped::isValid() const
{
	return _mappedView != NULL;
}


/// get file size
uint64_t MemoryMapped::size() const
{
	return _filesize;
}


/// get number of actually mapped bytes
size_t MemoryMapped::mappedSize() const
{
	return _mappedBytes;
}


/// replace mapping by a new one of the same file, offset MUST be a multiple of the page size
bool MemoryMapped::remap(uint64_t offset, size_t mappedBytes)
{
	_mappedView = mmap(nullptr, mappedBytes, PROT_WRITE | PROT_READ, MAP_SHARED, _file, 0);
	if (_mappedView == MAP_FAILED)
	{
		_mappedBytes = 0;
		_mappedView = NULL;
		return false;
	}

	_mappedBytes = mappedBytes;

	return true;
}


/// get OS page size (for remap)
int MemoryMapped::getpagesize()
{
#ifdef _WIN32
	SYSTEM_INFO sysInfo;
	GetSystemInfo(&sysInfo);
	return sysInfo.dwAllocationGranularity;
#else
	return sysconf(_SC_PAGESIZE); //::getpagesize();
#endif
}
