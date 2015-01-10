#include <gtest/gtest.h>
#include <Pixel.hpp>

TEST(PixelTest, operatorPixelAndPixel)
{
    ASSERT_EQ(Pixel(3,4,5), Pixel(1,1,1) + Pixel(2,3,4));
    ASSERT_EQ(Pixel(1,2,3), Pixel(2,3,4) - Pixel(1,1,1));
    ASSERT_EQ(Pixel(2,3,4), Pixel(1,1,1) * Pixel(2,3,4));
    ASSERT_EQ(Pixel(2,2,2), Pixel(4,6,8) / Pixel(2,3,4));
}

TEST(PixelTest, operatorPixelAndNumber)
{
    ASSERT_EQ(Pixel(3,3,3), Pixel(1,1,1) + 2);
    ASSERT_EQ(Pixel(0,1,2), Pixel(2,3,4) - 2);
    ASSERT_EQ(Pixel(3,3,3), Pixel(1,1,1) * 3);
    ASSERT_EQ(Pixel(2,3,4), Pixel(4,6,8) / 2);
}

TEST(PixelTest, lightness)
{
    ASSERT_EQ((4+6+8)/3, Pixel(4,6,8).Lightness());
}

TEST(PixelTest, range)
{
    ASSERT_EQ(Pixel(255,255,255), Pixel(0,0,0) - 1);
    ASSERT_EQ(Pixel(0,0,0), Pixel(255,255,255) + 1);
}

TEST(PixelTest, emptyPixel)
{
    Pixel pixel1, pixel2;
    ASSERT_TRUE(pixel1 == pixel2);
    pixel1.setRGB(0,0,0);
    ASSERT_FALSE(Pixel() == pixel1);
    pixel2.Label(1);
    ASSERT_FALSE(Pixel() == pixel2);

    ASSERT_FALSE(Pixel() == Pixel(0,0,0));
}