function FBRHS = FBLong(x, Ia, Ib, Ic, Id, Is, Ixi)
 global delta h sigma_X Ein tau0
 FBRHS = (0.10e2 / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * sqrt(0.2e1) * (x(5) ^ 16) + 0.60e2 * sqrt(0.2e1) * (sigma_X ^ 2) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 14) + 0.155e3 * sqrt(0.2e1) * (sigma_X ^ 4) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 12) + 0.225e3 * sqrt(0.2e1) * (sigma_X ^ 6) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 10) + 0.1605e4 / 0.8e1 * sqrt(0.2e1) * (sigma_X ^ 8) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 8) + 0.225e3 / 0.2e1 * sqrt(0.2e1) * (sigma_X ^ 10) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 6) + 0.155e3 / 0.4e1 * (sigma_X ^ 12) * sqrt(0.2e1) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 4) + 0.15e2 / 0.2e1 * (sigma_X ^ 14) * sqrt(0.2e1) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (x(5) ^ 2) + 0.5e1 / 0.8e1 * (sigma_X ^ 16) * sqrt(0.2e1) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4)) * x(1) ^ 2 + (0.192e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 + 0.128e3 * sigma_X * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) + 0.128e3 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * sigma_X * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.768e3 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 2) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2))) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 16) / 0.8e1 + (-0.128e3 * sigma_X * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) + 0.768e3 * (sigma_X ^ 3) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (tau0 ^ 2) + 0.256e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * x(6) * h ^ 2 - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (x(6) ^ 2) + 0.768e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 2) - 0.3072e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 2) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (x(6) ^ 2) + 0.960e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 2) + 0.3072e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 4) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.768e3 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * (sigma_X ^ 3) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) - 0.3072e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 2) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (tau0 ^ 2) - 0.768e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 2) * tau0 + 0.6144e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 2) * tau0) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 14) / 0.8e1 + (0.18432e5 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 4) * h * sqrt(pi) * x(4) * tau0 * x(6) + 0.1984e4 * (sigma_X ^ 5) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) + 0.1024e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * h * (sigma_X ^ 2) + 0.256e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * x(6) * h ^ 2 * (sigma_X ^ 2) - 0.512e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 2) + 0.1536e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 2) * tau0 + 0.1968e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 4) + 0.4992e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 6) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.3072e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 4) - 0.128e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (tau0 ^ 2) * (sigma_X ^ 2) - 0.1536e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 2) * (tau0 ^ 2) - 0.4096e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 2) * (tau0 ^ 3) - 0.4096e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 2) * tau0 - 0.9216e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 4) * h * sqrt(pi) * x(4) * (tau0 ^ 2) - 0.3072e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 4) * tau0 - 0.128e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 2) * (x(6) ^ 2) + 0.6144e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 2) * (tau0 ^ 2) - 0.9216e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 4) * h * sqrt(pi) * x(4) * (x(6) ^ 2) + 0.1984e4 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * (sigma_X ^ 5) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.512e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 2) * (tau0 ^ 3) - 0.768e3 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 3) + 0.1024e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 2) * (tau0 ^ 4)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 12) / 0.8e1 + (-0.8192e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 4) * tau0 + 0.21504e5 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 6) * h * sqrt(pi) * x(4) * tau0 * x(6) - 0.4608e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 4) * (tau0 ^ 2) + 0.4608e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 4) * tau0 + 0.12288e5 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 4) * (tau0 ^ 2) - 0.8192e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 4) * (tau0 ^ 3) - 0.1984e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 5) + 0.2880e4 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * (sigma_X ^ 7) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.2112e4 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 6) + 0.2880e4 * (sigma_X ^ 7) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 + 0.1056e4 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (tau0 ^ 2) * (sigma_X ^ 4) + 0.1056e4 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 4) * (x(6) ^ 2) - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 2) + 0.4224e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 8) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 2) * (tau0 ^ 4) - 0.768e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 2) * (x(6) ^ 2) * h ^ 2 * (sigma_X ^ 2) + 0.2048e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 4) * (tau0 ^ 4) - 0.4992e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 6) * tau0 + 0.4992e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 6) + 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * (x(6) ^ 3) * h ^ 2 * (sigma_X ^ 2) + 0.2048e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * h * (sigma_X ^ 4) + 0.1536e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 4) * (tau0 ^ 3) - 0.10752e5 * (sigma_X ^ 6) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (x(6) ^ 2) - 0.2112e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * x(6) * h ^ 2 * (sigma_X ^ 4) + 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 3) * x(6) * h ^ 2 * (sigma_X ^ 2) - 0.1536e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 4) - 0.10752e5 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 6) * h * sqrt(pi) * x(4) * (tau0 ^ 2)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 10) / 0.8e1 + (0.12288e5 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 8) * h * sqrt(pi) * x(4) * tau0 * x(6) - 0.5376e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 6) * (tau0 ^ 2) + 0.9216e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 6) * (tau0 ^ 2) - 0.6144e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 6) * (tau0 ^ 3) + 0.5376e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 6) * tau0 - 0.6144e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 6) * tau0 - 0.2880e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 7) + 0.2568e4 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * (sigma_X ^ 9) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.1248e4 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 8) + 0.2568e4 * (sigma_X ^ 9) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 - 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 4) * (tau0 ^ 4) - 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 4) + 0.1968e4 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 10) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.2976e4 * (sigma_X ^ 6) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (x(6) ^ 2) + 0.2976e4 * (sigma_X ^ 6) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (tau0 ^ 2) - 0.5952e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * x(6) * h ^ 2 * (sigma_X ^ 6) + 0.2048e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 3) * x(6) * h ^ 2 * (sigma_X ^ 4) - 0.3072e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 2) * (x(6) ^ 2) * h ^ 2 * (sigma_X ^ 4) - 0.6144e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 8) * h * sqrt(pi) * x(4) * (tau0 ^ 2) - 0.1792e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 6) + 0.1536e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 6) * (tau0 ^ 4) - 0.4224e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 8) * tau0 + 0.4224e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 8) + 0.2048e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * (x(6) ^ 3) * h ^ 2 * (sigma_X ^ 4) + 0.1536e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * h * (sigma_X ^ 6) + 0.1792e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 6) * (tau0 ^ 3) - 0.6144e4 * (sigma_X ^ 8) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (x(6) ^ 2)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 8) / 0.8e1 + (0.3072e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 8) * tau0 - 0.2048e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 8) * tau0 + 0.3072e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 8) * (tau0 ^ 2) + 0.3456e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 10) * h * sqrt(pi) * x(4) * tau0 * x(6) - 0.3072e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 8) * (tau0 ^ 2) - 0.2048e4 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 8) * (tau0 ^ 3) - 0.2568e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 9) + 0.1440e4 * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * delta * (sigma_X ^ 11) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.384e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 10) + 0.1440e4 * (sigma_X ^ 11) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 - 0.768e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 6) * (tau0 ^ 4) - 0.768e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 6) + 0.480e3 * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * (sigma_X ^ 12) * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.3264e4 * (sigma_X ^ 8) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (x(6) ^ 2) + 0.3264e4 * (sigma_X ^ 8) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (tau0 ^ 2) + 0.1024e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 8) * (tau0 ^ 3) - 0.1728e4 * (sigma_X ^ 10) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (x(6) ^ 2) - 0.1024e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 8) + 0.3072e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 3) * x(6) * h ^ 2 * (sigma_X ^ 6) - 0.6528e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * x(6) * h ^ 2 * (sigma_X ^ 8) - 0.4608e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 2) * (x(6) ^ 2) * h ^ 2 * (sigma_X ^ 6) - 0.1728e4 * (sigma_X ^ 10) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (tau0 ^ 2) + 0.512e3 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * h * (sigma_X ^ 8) + 0.512e3 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 8) * (tau0 ^ 4) + 0.1968e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 10) + 0.3072e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * (x(6) ^ 3) * h ^ 2 * (sigma_X ^ 6) - 0.1968e4 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 10) * tau0) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 6) / 0.8e1 + (0.864e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 10) * tau0 - 0.256e3 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 10) * (tau0 ^ 3) + 0.384e3 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 10) * (tau0 ^ 2) + 0.384e3 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 12) * h * sqrt(pi) * x(4) * tau0 * x(6) - 0.864e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 10) * (tau0 ^ 2) - 0.256e3 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 10) * tau0 - 0.1440e4 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 11) + 0.496e3 * (sigma_X ^ 13) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 + 0.48e2 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 12) + 0.496e3 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * delta * (sigma_X ^ 13) - 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 8) * (tau0 ^ 4) + 0.1728e4 * (sigma_X ^ 10) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (x(6) ^ 2) - 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 8) + 0.48e2 * (sigma_X ^ 14) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(4) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) + 0.1728e4 * (sigma_X ^ 10) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (tau0 ^ 2) + 0.64e2 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * h * (sigma_X ^ 10) + 0.64e2 * sqrt(pi) * x(4) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 10) * (tau0 ^ 4) - 0.3072e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 2) * (x(6) ^ 2) * h ^ 2 * (sigma_X ^ 8) + 0.2048e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * (x(6) ^ 3) * h ^ 2 * (sigma_X ^ 8) + 0.2048e4 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 3) * x(6) * h ^ 2 * (sigma_X ^ 8) - 0.480e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 12) * tau0 + 0.288e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 10) * (tau0 ^ 3) - 0.192e3 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (sigma_X ^ 12) * h * sqrt(pi) * x(4) * (x(6) ^ 2) - 0.288e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 10) - 0.3456e4 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * tau0 * (sigma_X ^ 10) * x(6) + 0.480e3 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 12) - 0.192e3 * (sigma_X ^ 12) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(4) * (tau0 ^ 2)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 4) / 0.8e1 + (0.96e2 * (sigma_X ^ 15) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 + 0.32e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 12) * (tau0 ^ 3) + 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 3) * x(6) * h ^ 2 * (sigma_X ^ 10) + 0.416e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (tau0 ^ 2) * (sigma_X ^ 12) + 0.96e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 2) * h * (sigma_X ^ 12) * tau0 + 0.48e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 14) + 0.96e2 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * delta * (sigma_X ^ 15) - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 4) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 10) - 0.32e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * (x(6) ^ 3) * h * (sigma_X ^ 12) - 0.768e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * (tau0 ^ 2) * (x(6) ^ 2) * h ^ 2 * (sigma_X ^ 10) - 0.96e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * x(6) * h * (sigma_X ^ 12) * (tau0 ^ 2) - 0.48e2 * sqrt(pi) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) * exp(-((-tau0 + x(6)) ^ 2 / (x(5) ^ 2 + sigma_X ^ 2))) * h * (sigma_X ^ 14) * tau0 - 0.128e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * h ^ 2 * (sigma_X ^ 10) * (tau0 ^ 4) + 0.416e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 12) * (x(6) ^ 2) - 0.496e3 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 13) - 0.832e3 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * tau0 * (sigma_X ^ 12) * x(6) + 0.512e3 * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * tau0 * (x(6) ^ 3) * h ^ 2 * (sigma_X ^ 10)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X * (x(5) ^ 2) / 0.8e1 + (0.8e1 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * delta * (sigma_X ^ 17) + 0.32e2 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (tau0 ^ 2) * (sigma_X ^ 14) - 0.64e2 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * tau0 * (sigma_X ^ 14) * x(6) + 0.8e1 * (sigma_X ^ 17) * sqrt(pi) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * x(3) ^ 2 + 0.32e2 * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * exp(-(2 * (-tau0 + x(6)) ^ 2 / (2 * x(5) ^ 2 + sigma_X ^ 2))) * h ^ 2 * (sigma_X ^ 14) * (x(6) ^ 2) - 0.96e2 * sqrt((2 * x(5) ^ 2 + sigma_X ^ 2)) * sqrt((x(5) ^ 2 + sigma_X ^ 2)) * sqrt(pi) * (sigma_X ^ 15)) * pi ^ (-0.1e1 / 0.2e1) * ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) * ((x(5) ^ 2 + sigma_X ^ 2) ^ (-0.9e1 / 0.2e1)) / sigma_X / 0.8e1 - (1 / (2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4 / (x(5) ^ 2 + sigma_X ^ 2) ^ 4 * sigma_X ^ 16 / x(5) ^ 2) + (-0.12e2 / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * Ia * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 15) - 0.72e2 * Ia * (sigma_X ^ 2) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 13) - 0.186e3 * Ia * (sigma_X ^ 4) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 11) - 0.270e3 * Ia * (sigma_X ^ 6) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 9) - 0.963e3 / 0.4e1 * Ia * (sigma_X ^ 8) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 7) - 0.135e3 * Ia * (sigma_X ^ 10) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 5) - 0.93e2 / 0.2e1 * (sigma_X ^ 12) * Ia / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 3) - 0.9e1 * (sigma_X ^ 14) * Ia / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * x(5) - 0.3e1 / 0.4e1 / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * Ia * (sigma_X ^ 16) * pi ^ (-0.1e1 / 0.2e1) / x(5)) / x(1) + (0.8e1 * Is / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 16) + 0.16e2 / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 15) + 0.48e2 * Is * (sigma_X ^ 2) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 14) + 0.96e2 * (sigma_X ^ 2) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 13) + 0.124e3 * Is * (sigma_X ^ 4) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 12) + 0.248e3 * (sigma_X ^ 4) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 11) + 0.180e3 * Is * (sigma_X ^ 6) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 10) + 0.360e3 * (sigma_X ^ 6) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 9) + 0.321e3 / 0.2e1 * Is * (sigma_X ^ 8) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 8) + 0.321e3 * (sigma_X ^ 8) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 7) + 0.90e2 * Is * (sigma_X ^ 10) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 6) + 0.180e3 * (sigma_X ^ 10) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 5) + 0.31e2 * Is * (sigma_X ^ 12) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 4) + 0.62e2 * (sigma_X ^ 12) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 3) + 0.6e1 * Is * (sigma_X ^ 14) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * (x(5) ^ 2) + 0.12e2 * (sigma_X ^ 14) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) * x(3) * Ic / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) * x(5) + Is * (sigma_X ^ 16) / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * pi ^ (-0.1e1 / 0.2e1) / 0.2e1 + 0.1e1 / ((2 * x(5) ^ 2 + sigma_X ^ 2) ^ 4) / ((x(5) ^ 2 + sigma_X ^ 2) ^ 4) * (sigma_X ^ 16) * x(3) * Ic * pi ^ (-0.1e1 / 0.2e1) / x(5)) / x(1) ^ 2;

end
