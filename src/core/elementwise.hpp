// Extend std::vector<double> with elementwise operations.
// This is useful for fast data manipulation
//
// ------------------------------------------------------------------------------
// Author:       Daniel Winney (2022)
// Affiliation:  Joint Physics Analysis Center (JPAC),
//               South China Normal Univeristy (SCNU)
// Email:        dwinney@iu.alumni.edu
// ------------------------------------------------------------------------------

#ifndef ELEMENTWISE_HPP
#define ELEMENTWISE_HPP

namespace jpacPhoto
{
    
    // ---------------------------------------------------------------------------
    // Element-wise operations on data vectors

    inline std::vector<double>operator+(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] + rhs[i] );
        };
        return result;
    };

    inline std::vector<double>operator-(std::vector<double> lhs, std::vector<double> rhs)
    {
        if (lhs.size() != rhs.size()) return error("Attempted to add two vectors of different sizes!", std::vector<double>());

        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i] - rhs[i] );
        };
        return result;
    };

    // Given two vector<double>s of the same size, calculate the average element wise
    inline std::vector<double>operator*( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]*c );
        };
        return result;
    };

    inline std::vector<double>operator*(double c, std::vector<double> rhs)
    {
        std::vector<double> result;
        for (int i = 0; i < rhs.size(); i++)
        {
            result.push_back( c*rhs[i] );
        };
        return result;
    };

    inline std::vector<double>operator/( std::vector<double> lhs, double c)
    {
        std::vector<double> result;
        for (int i = 0; i < lhs.size(); i++)
        {
            result.push_back( lhs[i]/c );
        };
        return result;
    };
};

#endif