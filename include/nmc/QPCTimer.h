// QPCTimer.h
// QueryPerformanceFrequency と QueryPerformanceCounter
// を利用したタイマークラス
#ifdef _WIN32

#pragma once
#include <windows.h>

class QPCTimer{
public:
  QPCTimer(){
    LARGE_INTEGER li;
    BOOL re = QueryPerformanceFrequency(&li);
    assert(re);
    freq_ = double(li.QuadPart / LONGLONG(1000));
    start_();
  }

  ~QPCTimer(){};

  void restart(){
    start_();
  }

  double elapsed() const{
    LARGE_INTEGER li;
    QueryPerformanceCounter(&li);
    return double(li.QuadPart - begin_) / freq_;
  }

private:
  void start_(){
    LARGE_INTEGER li;
    BOOL re = QueryPerformanceCounter(&li);
    assert(re);
    begin_ = li.QuadPart;
  }

private:
  double   freq_;
  LONGLONG begin_;
};

#endif //_WIN32