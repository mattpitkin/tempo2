
namespace testing {
namespace internal {
enum GTestColor {
  COLOR_DEFAULT,
  COLOR_RED,
  COLOR_GREEN,
  COLOR_YELLOW,
  COLOR_CYAN
};
void ColoredPrintf(GTestColor color, const char* fmt, ...);
}
}


#define NOTE(...) testing::internal::ColoredPrintf(testing::internal::COLOR_CYAN,"[ NOTE     ] "); printf(__VA_ARGS__)
