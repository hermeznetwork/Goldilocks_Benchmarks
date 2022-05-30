# Build Instructions
## Pre-requisites 
### build-essential, cmake and OMP
```
$ apt install -y build-essential libomp-dev cmake
```
### Install google benchmark

```
$ git clone https://github.com/google/benchmark.git
$ cd benchmark
$ cmake -E make_directory "build"
$ cmake -E chdir "build" cmake -DBENCHMARK_DOWNLOAD_DEPENDENCIES=on -DCMAKE_BUILD_TYPE=Release ../
$ cmake --build "build" --config Release
```

Run

```
& ./build/benchmark --benchmark_counters_tabular=true
```

## Results
### AMD Ryzen Threadripper 3990X 64-Core Processor & 256GB RAM

```bash
022-05-30T11:03:26+02:00
Running ./build/benchmark
Run on (128 X 3599.63 MHz CPU s)
CPU Caches:
  L1 Data 32 KiB (x64)
  L1 Instruction 32 KiB (x64)
  L2 Unified 512 KiB (x64)
  L3 Unified 16384 KiB (x16)
Load Average: 9.54, 12.67, 13.50
-------------------------------------------------------------------------------------------------------
Benchmark                                   Time             CPU   Iterations BytesProcessed       Rate
-------------------------------------------------------------------------------------------------------
POSEIDON_BENCH/1/real_time              29990 us        29992 us           23     30.5275M/s  2.99903us
POSEIDON_BENCH/2/real_time              15087 us        15088 us           46     60.6819M/s  3.01747us
POSEIDON_BENCH/4/real_time               7603 us         7603 us           93     120.424M/s  3.04102us
POSEIDON_BENCH/8/real_time               3875 us         3875 us          182     236.263M/s  3.10003us
POSEIDON_BENCH/16/real_time              1950 us         1950 us          360     469.436M/s  3.12043us
POSEIDON_BENCH/32/real_time               994 us          994 us          533     921.016M/s  3.18093us
POSEIDON_BENCH/64/real_time              1052 us         1053 us          630     869.874M/s  6.73589us
POSEIDON_BENCH/128/real_time              929 us          842 us          554     985.358M/s  11.8929us
POSEIDON_BENCH/60/real_time              1003 us         1003 us          590     912.982M/s  6.01673us
POSEIDON_BENCH/62/real_time               581 us          581 us          928     1.53998G/s  3.59956us
POSEIDON_BENCH/64/real_time               967 us          967 us          609      946.88M/s  6.18809us
POSEIDON_BENCH/66/real_time              1027 us         1028 us          692     891.064M/s   6.7812us
POSEIDON_BENCH/68/real_time              1010 us         1010 us          708     906.583M/s  6.86709us
LINEAR_HASH_SINGLE_BENCH/real_time       39.7 us         39.7 us        16696     19.2061M/s  3.05567us
LINEAR_HASH_BENCH/60/real_time           25.7 s          25.7 s             1     996.084M/s   3.5351us
LINEAR_HASH_BENCH/62/real_time           25.4 s          25.4 s             1      1007.1M/s  3.61299us
LINEAR_HASH_BENCH/64/real_time           24.7 s          24.7 s             1     1036.17M/s   3.6249us
LINEAR_HASH_BENCH/66/real_time           33.8 s          24.3 s             1     756.553M/s  5.11978us
LINEAR_HASH_BENCH/68/real_time           33.0 s          23.6 s             1     776.542M/s  5.13914us
LINEAR_HASH_BENCH/128/real_time          25.8 s          25.6 s             1     991.897M/s  7.57339us
MERKLE_TREE_BENCH/60/real_time           28.3 s          28.2 s             1     903.509M/s  3.89731us
MERKLE_TREE_BENCH/62/real_time           28.3 s          28.3 s             1     905.153M/s  4.01991us
MERKLE_TREE_BENCH/64/real_time           34.3 s          27.9 s             1     746.731M/s  5.02993us
MERKLE_TREE_BENCH/66/real_time           36.0 s          26.6 s             1     710.177M/s  5.45412us
MERKLE_TREE_BENCH/68/real_time           35.2 s          25.8 s             1     727.157M/s  5.48817us
MERKLE_TREE_BENCH/128/real_time          28.4 s          28.4 s             1      902.95M/s  8.31942us
iNTT_BENCH/60/real_time                  4.54 s          4.11 s             1     1.37592G/s   2.72545s
iNTT_BENCH/62/real_time                  4.70 s          4.46 s             1     1.32846G/s    2.9169s
iNTT_BENCH/64/real_time                  4.60 s          4.37 s             1     1.35977G/s   2.94168s
iNTT_BENCH/66/real_time                  5.00 s          4.58 s             1     1.25013G/s   3.29965s
iNTT_BENCH/68/real_time                  4.80 s          4.50 s             1     1.30182G/s   3.26466s
iNTT_BENCH/128/real_time                 3.47 s          3.44 s             1     1.80228G/s   4.43883s
-----------------------------------------------------------------------------
Benchmark                                   Time             CPU   Iterations
-----------------------------------------------------------------------------
LDE_BENCH/60/real_time                   14.7 s          14.3 s             1
LDE_BENCH/62/real_time                   14.3 s          13.5 s             1
LDE_BENCH/64/real_time                   14.4 s          13.9 s             1
LDE_BENCH/66/real_time                   14.1 s          12.6 s             1
LDE_BENCH/68/real_time                   14.3 s          12.7 s             1
LDE_BENCH/128/real_time                  11.2 s          10.5 s             1
```


## Linear HASH
```
---------
(0,0): 0x1              (0,1): 0x2              (0,2): 0x3              (0,3): 0x4              (0,4): 0x5              (0,5): 0x6              (0,6): 0x7              (0,7): 0x8              (0,8): 0x9              (0,9): 0xa
(1,0): 0x1              (1,1): 0x2              (1,2): 0x3              (1,3): 0x4              (1,4): 0x5              (1,5): 0x6              (1,6): 0x7              (1,7): 0x8              (1,8): 0x9              (1,9): 0xa
(2,0): 0x2              (2,1): 0x4              (2,2): 0x6              (2,3): 0x8              (2,4): 0xa              (2,5): 0xc              (2,6): 0xe              (2,7): 0x10             (2,8): 0x12             (2,9): 0x14
(3,0): 0x3              (3,1): 0x6              (3,2): 0x9              (3,3): 0xc              (3,4): 0xf              (3,5): 0x12             (3,6): 0x15             (3,7): 0x18             (3,8): 0x1b             (3,9): 0x1e
(4,0): 0x5              (4,1): 0xa              (4,2): 0xf              (4,3): 0x14             (4,4): 0x19             (4,5): 0x1e             (4,6): 0x23             (4,7): 0x28             (4,8): 0x2d             (4,9): 0x32
(5,0): 0x8              (5,1): 0x10             (5,2): 0x18             (5,3): 0x20             (5,4): 0x28             (5,5): 0x30             (5,6): 0x38             (5,7): 0x40             (5,8): 0x48             (5,9): 0x50
(6,0): 0xd              (6,1): 0x1a             (6,2): 0x27             (6,3): 0x34             (6,4): 0x41             (6,5): 0x4e             (6,6): 0x5b             (6,7): 0x68             (6,8): 0x75             (6,9): 0x82
(7,0): 0x15             (7,1): 0x2a             (7,2): 0x3f             (7,3): 0x54             (7,4): 0x69             (7,5): 0x7e             (7,6): 0x93             (7,7): 0xa8             (7,8): 0xbd             (7,9): 0xd2
(8,0): 0x22             (8,1): 0x44             (8,2): 0x66             (8,3): 0x88             (8,4): 0xaa             (8,5): 0xcc             (8,6): 0xee             (8,7): 0x110            (8,8): 0x132            (8,9): 0x154
(9,0): 0x37             (9,1): 0x6e             (9,2): 0xa5             (9,3): 0xdc             (9,4): 0x113            (9,5): 0x14a            (9,6): 0x181            (9,7): 0x1b8            (9,8): 0x1ef            (9,9): 0x226
(10,0): 0x59            (10,1): 0xb2            (10,2): 0x10b           (10,3): 0x164           (10,4): 0x1bd           (10,5): 0x216           (10,6): 0x26f           (10,7): 0x2c8           (10,8): 0x321           (10,9): 0x37a
(11,0): 0x90            (11,1): 0x120           (11,2): 0x1b0           (11,3): 0x240           (11,4): 0x2d0           (11,5): 0x360           (11,6): 0x3f0           (11,7): 0x480           (11,8): 0x510           (11,9): 0x5a0
(12,0): 0xe9            (12,1): 0x1d2           (12,2): 0x2bb           (12,3): 0x3a4           (12,4): 0x48d           (12,5): 0x576           (12,6): 0x65f           (12,7): 0x748           (12,8): 0x831           (12,9): 0x91a
(13,0): 0x179           (13,1): 0x2f2           (13,2): 0x46b           (13,3): 0x5e4           (13,4): 0x75d           (13,5): 0x8d6           (13,6): 0xa4f           (13,7): 0xbc8           (13,8): 0xd41           (13,9): 0xeba
(14,0): 0x262           (14,1): 0x4c4           (14,2): 0x726           (14,3): 0x988           (14,4): 0xbea           (14,5): 0xe4c           (14,6): 0x10ae          (14,7): 0x1310          (14,8): 0x1572          (14,9): 0x17d4
(15,0): 0x3db           (15,1): 0x7b6           (15,2): 0xb91           (15,3): 0xf6c           (15,4): 0x1347          (15,5): 0x1722          (15,6): 0x1afd          (15,7): 0x1ed8          (15,8): 0x22b3          (15,9): 0x268e
---------
linear_hash ( 0x1  0x2  0x3  0x4  0x5  0x6  0x7  0x8  0  0  0  0 ) -> (  0xd110aa6a46373941  0x8f238fcceb658894  0x9cd4f8353866fb4f  0x274913f0007aa232 )
linear_hash ( 0x9  0xa  0  0  0  0  0  0  0xd110aa6a46373941  0x8f238fcceb658894  0x9cd4f8353866fb4f  0x274913f0007aa232 ) -> (  0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9 )
result:  0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9 
---------
linear_hash ( 0x1  0x2  0x3  0x4  0x5  0x6  0x7  0x8  0  0  0  0 ) -> (  0xd110aa6a46373941  0x8f238fcceb658894  0x9cd4f8353866fb4f  0x274913f0007aa232 )
linear_hash ( 0x9  0xa  0  0  0  0  0  0  0xd110aa6a46373941  0x8f238fcceb658894  0x9cd4f8353866fb4f  0x274913f0007aa232 ) -> (  0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9 )
result:  0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9 
---------
linear_hash ( 0x2  0x4  0x6  0x8  0xa  0xc  0xe  0x10  0  0  0  0 ) -> (  0xb9286454615de76c  0x618c71e1780f9e40  0xd0c3d45afbef68e  0x817238a7d7375823 )
linear_hash ( 0x12  0x14  0  0  0  0  0  0  0xb9286454615de76c  0x618c71e1780f9e40  0xd0c3d45afbef68e  0x817238a7d7375823 ) -> (  0x887d3121ac84d143  0xb24e03a8d54216e9  0xdd8ffde523c87f63  0x125730aa2ec88be8 )
result:  0x887d3121ac84d143  0xb24e03a8d54216e9  0xdd8ffde523c87f63  0x125730aa2ec88be8 
---------
linear_hash ( 0x3  0x6  0x9  0xc  0xf  0x12  0x15  0x18  0  0  0  0 ) -> (  0xfe61611807751319  0x5b97cd7b2f836067  0xd32b09819862c264  0xae4f2ae427267c92 )
linear_hash ( 0x1b  0x1e  0  0  0  0  0  0  0xfe61611807751319  0x5b97cd7b2f836067  0xd32b09819862c264  0xae4f2ae427267c92 ) -> (  0xc0a7b2dec34af7d4  0x710202c8301bbc66  0x6770df66b764e649  0xebdbb06c751fb842 )
result:  0xc0a7b2dec34af7d4  0x710202c8301bbc66  0x6770df66b764e649  0xebdbb06c751fb842 
---------
linear_hash ( 0x5  0xa  0xf  0x14  0x19  0x1e  0x23  0x28  0  0  0  0 ) -> (  0xd62a43b76be50633  0xf535ee0a4518b351  0x67f60264cde11ffa  0x64a917418878d2e9 )
linear_hash ( 0x2d  0x32  0  0  0  0  0  0  0xd62a43b76be50633  0xf535ee0a4518b351  0x67f60264cde11ffa  0x64a917418878d2e9 ) -> (  0x5c44254efc9c737c  0x8342f0df861b8f42  0xa69eb5140f15574d  0x35223f41d305436f )
result:  0x5c44254efc9c737c  0x8342f0df861b8f42  0xa69eb5140f15574d  0x35223f41d305436f 
---------
linear_hash ( 0x8  0x10  0x18  0x20  0x28  0x30  0x38  0x40  0  0  0  0 ) -> (  0xd50f30d0097c9fde  0xf5afb5e0d8dac1eb  0x5e4cb1a9972fc6dd  0x7ff8090d73481c86 )
linear_hash ( 0x48  0x50  0  0  0  0  0  0  0xd50f30d0097c9fde  0xf5afb5e0d8dac1eb  0x5e4cb1a9972fc6dd  0x7ff8090d73481c86 ) -> (  0x97d55f1b6e238d9e  0xbeddd5f9a1319adf  0xba2c3c4f4072d58b  0xd1dc1e4a71e14750 )
result:  0x97d55f1b6e238d9e  0xbeddd5f9a1319adf  0xba2c3c4f4072d58b  0xd1dc1e4a71e14750 
---------
linear_hash ( 0xd  0x1a  0x27  0x34  0x41  0x4e  0x5b  0x68  0  0  0  0 ) -> (  0x9b8e394e594911fd  0x2b4f622fc27da3e4  0x2bb705ab95b19123  0xf88b3e8a4d786885 )
linear_hash ( 0x75  0x82  0  0  0  0  0  0  0x9b8e394e594911fd  0x2b4f622fc27da3e4  0x2bb705ab95b19123  0xf88b3e8a4d786885 ) -> (  0xea6c029509ea209d  0x6023440929b8a71b  0x240d3963ea04cb3a  0x898007e5765208a4 )
result:  0xea6c029509ea209d  0x6023440929b8a71b  0x240d3963ea04cb3a  0x898007e5765208a4 
---------
linear_hash ( 0x15  0x2a  0x3f  0x54  0x69  0x7e  0x93  0xa8  0  0  0  0 ) -> (  0x19ca61218e68bf9f  0xc7db13d8c3e44c10  0xf3cb9f0dd2932a4b  0x5ace7010f0e31bed )
linear_hash ( 0xbd  0xd2  0  0  0  0  0  0  0x19ca61218e68bf9f  0xc7db13d8c3e44c10  0xf3cb9f0dd2932a4b  0x5ace7010f0e31bed ) -> (  0xfa9e9de778ca8f2b  0x178c1550469d3e04  0x560a99f234daf023  0xceb5c4f7a4828ada )
result:  0xfa9e9de778ca8f2b  0x178c1550469d3e04  0x560a99f234daf023  0xceb5c4f7a4828ada 
---------
linear_hash ( 0x22  0x44  0x66  0x88  0xaa  0xcc  0xee  0x110  0  0  0  0 ) -> (  0x9492a13254db601c  0x44f475695f5f1f4d  0xb3fe7c17afff7353  0x552024eb2ed4d89d )
linear_hash ( 0x132  0x154  0  0  0  0  0  0  0x9492a13254db601c  0x44f475695f5f1f4d  0xb3fe7c17afff7353  0x552024eb2ed4d89d ) -> (  0xb9d0b9f50b3f55b  0x3fb21df04c1f5627  0x155d9ab00534358c  0xd6bd9266e7fe968e )
result:  0xb9d0b9f50b3f55b  0x3fb21df04c1f5627  0x155d9ab00534358c  0xd6bd9266e7fe968e 
---------
linear_hash ( 0x37  0x6e  0xa5  0xdc  0x113  0x14a  0x181  0x1b8  0  0  0  0 ) -> (  0x5bb96798686c8903  0x3cdeb4f620463dad  0x8e5bb61155675a3e  0xfc6be8da247d43b1 )
linear_hash ( 0x1ef  0x226  0  0  0  0  0  0  0x5bb96798686c8903  0x3cdeb4f620463dad  0x8e5bb61155675a3e  0xfc6be8da247d43b1 ) -> (  0x3b0d55bb5ca697f4  0x94ef98c6a39b6734  0x26f595f0ef381db3  0x5fd1278619025f2a )
result:  0x3b0d55bb5ca697f4  0x94ef98c6a39b6734  0x26f595f0ef381db3  0x5fd1278619025f2a 
---------
linear_hash ( 0x59  0xb2  0x10b  0x164  0x1bd  0x216  0x26f  0x2c8  0  0  0  0 ) -> (  0xd8a5f3a7bdc02976  0x5534eddc8ccb6eaa  0xc17816c0c4a6a980  0xc74fc43832a053ea )
linear_hash ( 0x321  0x37a  0  0  0  0  0  0  0xd8a5f3a7bdc02976  0x5534eddc8ccb6eaa  0xc17816c0c4a6a980  0xc74fc43832a053ea ) -> (  0xe424a2fc736a1a65  0xc621e67a0b7b00b7  0x473a7ef505c7d1c6  0x9df81ccc4e06d82e )
result:  0xe424a2fc736a1a65  0xc621e67a0b7b00b7  0x473a7ef505c7d1c6  0x9df81ccc4e06d82e 
---------
linear_hash ( 0x90  0x120  0x1b0  0x240  0x2d0  0x360  0x3f0  0x480  0  0  0  0 ) -> (  0xf35dcbb276ca7871  0x29b8157c62b278ef  0x3802e8a79b9e3be6  0xff25ef390c09b3d4 )
linear_hash ( 0x510  0x5a0  0  0  0  0  0  0  0xf35dcbb276ca7871  0x29b8157c62b278ef  0x3802e8a79b9e3be6  0xff25ef390c09b3d4 ) -> (  0xe3b9b669a9211fd7  0xce9169262798cf4b  0xac9969f47d1124bf  0x1bffc4d7ea615df4 )
result:  0xe3b9b669a9211fd7  0xce9169262798cf4b  0xac9969f47d1124bf  0x1bffc4d7ea615df4 
---------
linear_hash ( 0xe9  0x1d2  0x2bb  0x3a4  0x48d  0x576  0x65f  0x748  0  0  0  0 ) -> (  0x6ef4d570af8689f2  0x673f144c14cfc80a  0x300d4914408e37a6  0x65463a817c871952 )
linear_hash ( 0x831  0x91a  0  0  0  0  0  0  0x6ef4d570af8689f2  0x673f144c14cfc80a  0x300d4914408e37a6  0x65463a817c871952 ) -> (  0x4c722ebc0c780dc1  0x9387e7b8d8709838  0x2d25676159e7a977  0x872e8e299f308877 )
result:  0x4c722ebc0c780dc1  0x9387e7b8d8709838  0x2d25676159e7a977  0x872e8e299f308877 
---------
linear_hash ( 0x179  0x2f2  0x46b  0x5e4  0x75d  0x8d6  0xa4f  0xbc8  0  0  0  0 ) -> (  0x1ba9b2dc92396d55  0x2b5199ecb089aa  0xdae1be9bacaef872  0x7392280c72985e39 )
linear_hash ( 0xd41  0xeba  0  0  0  0  0  0  0x1ba9b2dc92396d55  0x2b5199ecb089aa  0xdae1be9bacaef872  0x7392280c72985e39 ) -> (  0xbc3b546f83d7109b  0x644423ca0a26db0b  0x1ebf1e2b1fceec1e  0x37de800ea639f821 )
result:  0xbc3b546f83d7109b  0x644423ca0a26db0b  0x1ebf1e2b1fceec1e  0x37de800ea639f821 
---------
linear_hash ( 0x262  0x4c4  0x726  0x988  0xbea  0xe4c  0x10ae  0x1310  0  0  0  0 ) -> (  0xda0b9b1dd4d11056  0x6d68edc65698a6c4  0x2dd5d0474dbab001  0xaea9d8f4086f4139 )
linear_hash ( 0x1572  0x17d4  0  0  0  0  0  0  0xda0b9b1dd4d11056  0x6d68edc65698a6c4  0x2dd5d0474dbab001  0xaea9d8f4086f4139 ) -> (  0x459a0c9d2b83f3fd  0xa058e4a417f8daa9  0xf09d63da10517ee1  0x9a1ca1937d2d018f )
result:  0x459a0c9d2b83f3fd  0xa058e4a417f8daa9  0xf09d63da10517ee1  0x9a1ca1937d2d018f 
---------
linear_hash ( 0x3db  0x7b6  0xb91  0xf6c  0x1347  0x1722  0x1afd  0x1ed8  0  0  0  0 ) -> (  0xccfa35c3ed7cab0c  0x6832774945a56705  0xae14b5c33b1043f7  0xff8a6fd855c65178 )
linear_hash ( 0x22b3  0x268e  0  0  0  0  0  0  0xccfa35c3ed7cab0c  0x6832774945a56705  0xae14b5c33b1043f7  0xff8a6fd855c65178 ) -> (  0x53e60202c1f7379d  0x7537f186259879ba  0xa491aedb479d0e7b  0xb1bc4f02ac188223 )
result:  0x53e60202c1f7379d  0x7537f186259879ba  0xa491aedb479d0e7b  0xb1bc4f02ac188223 
---------
```

```
Merkle Tree
##########
hash ( 0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9  0x52a855f38bb05faf  0xa31e6272e294d7ab  0xd5499446c76d4674  0xaec3c33bd1409ae9  0  0  0  0 ) -> 
        ( 0x3f4b5ee8ca6e3a5e  0x473dacc190771f62  0x75f5d02e268d2281  0xfda76ae0848640e1 )
hash ( 0x887d3121ac84d143  0xb24e03a8d54216e9  0xdd8ffde523c87f63  0x125730aa2ec88be8  0xc0a7b2dec34af7d4  0x710202c8301bbc66  0x6770df66b764e649  0xebdbb06c751fb842  0  0  0  0 ) -> 
        ( 0xcd29fda5f6215c51  0x66efe8683d9d1c46  0x1fc267ba156f9b8d  0x7a6afe36d005b503 )
hash ( 0x5c44254efc9c737c  0x8342f0df861b8f42  0xa69eb5140f15574d  0x35223f41d305436f  0x97d55f1b6e238d9e  0xbeddd5f9a1319adf  0xba2c3c4f4072d58b  0xd1dc1e4a71e14750  0  0  0  0 ) -> 
        ( 0xdd2b1ebd6a83a3e6  0x59fc5c6b36f37d3  0x60037c5529c90946  0xdd0555eb45b555cf )
hash ( 0xea6c029509ea209d  0x6023440929b8a71b  0x240d3963ea04cb3a  0x898007e5765208a4  0xfa9e9de778ca8f2b  0x178c1550469d3e04  0x560a99f234daf023  0xceb5c4f7a4828ada  0  0  0  0 ) -> 
        ( 0xf317d3ab2a2dde6e  0x6aecf4a54c5f03f1  0xdf08fce16f634ecd  0x3723e8307dd223d4 )
hash ( 0xb9d0b9f50b3f55b  0x3fb21df04c1f5627  0x155d9ab00534358c  0xd6bd9266e7fe968e  0x3b0d55bb5ca697f4  0x94ef98c6a39b6734  0x26f595f0ef381db3  0x5fd1278619025f2a  0  0  0  0 ) -> 
        ( 0x46428d3a7e3af2ad  0xfcb12af8a292231  0xeaf048819ddd6685  0x8a09b9f2e7b46655 )
hash ( 0xe424a2fc736a1a65  0xc621e67a0b7b00b7  0x473a7ef505c7d1c6  0x9df81ccc4e06d82e  0xe3b9b669a9211fd7  0xce9169262798cf4b  0xac9969f47d1124bf  0x1bffc4d7ea615df4  0  0  0  0 ) -> 
        ( 0xdbe6a4ea2aed40a6  0xa0592f77881ce6aa  0xc65a926165883e28  0x22e17bfe0ddb1619 )
hash ( 0x4c722ebc0c780dc1  0x9387e7b8d8709838  0x2d25676159e7a977  0x872e8e299f308877  0xbc3b546f83d7109b  0x644423ca0a26db0b  0x1ebf1e2b1fceec1e  0x37de800ea639f821  0  0  0  0 ) -> 
        ( 0x325b20fae3dcd14b  0x9e82fac8872f46a  0x2195fbc523db2e7  0xe9d58cb4cc97632a )
hash ( 0x459a0c9d2b83f3fd  0xa058e4a417f8daa9  0xf09d63da10517ee1  0x9a1ca1937d2d018f  0x53e60202c1f7379d  0x7537f186259879ba  0xa491aedb479d0e7b  0xb1bc4f02ac188223  0  0  0  0 ) -> 
        ( 0xf4a0e645b580c67f  0x4d94af202c6828ce  0x424a5dfe9d958fc9  0xe8b63d2be82a79cc )
##########
hash ( 0x3f4b5ee8ca6e3a5e  0x473dacc190771f62  0x75f5d02e268d2281  0xfda76ae0848640e1  0xcd29fda5f6215c51  0x66efe8683d9d1c46  0x1fc267ba156f9b8d  0x7a6afe36d005b503  0  0  0  0 ) -> 
        ( 0xcd4e5de2ec9bc6b4  0xf888a029405b8f78  0xe780a71cf5cffdc2  0x6242ff456defee34 )
hash ( 0xdd2b1ebd6a83a3e6  0x59fc5c6b36f37d3  0x60037c5529c90946  0xdd0555eb45b555cf  0xf317d3ab2a2dde6e  0x6aecf4a54c5f03f1  0xdf08fce16f634ecd  0x3723e8307dd223d4  0  0  0  0 ) -> 
        ( 0x1bc59446b9ed8a96  0x4ecb0c09df8494ca  0x8cdedaa542e11cd8  0x8180763c48eff0c6 )
hash ( 0x46428d3a7e3af2ad  0xfcb12af8a292231  0xeaf048819ddd6685  0x8a09b9f2e7b46655  0xdbe6a4ea2aed40a6  0xa0592f77881ce6aa  0xc65a926165883e28  0x22e17bfe0ddb1619  0  0  0  0 ) -> 
        ( 0x6f159dfd0835443d  0x55d8a9877b672439  0x327dc934599e4bb7  0xa50a8c9661151fe5 )
hash ( 0x325b20fae3dcd14b  0x9e82fac8872f46a  0x2195fbc523db2e7  0xe9d58cb4cc97632a  0xf4a0e645b580c67f  0x4d94af202c6828ce  0x424a5dfe9d958fc9  0xe8b63d2be82a79cc  0  0  0  0 ) -> 
        ( 0x3b1224409e600901  0x2f2e2d34681eb87b  0x28de89ca8e928bc2  0xd3ec419ff5954d92 )
##########
hash ( 0xcd4e5de2ec9bc6b4  0xf888a029405b8f78  0xe780a71cf5cffdc2  0x6242ff456defee34  0x1bc59446b9ed8a96  0x4ecb0c09df8494ca  0x8cdedaa542e11cd8  0x8180763c48eff0c6  0  0  0  0 ) -> 
        ( 0x2e6407ab14535ea  0xa07ff473834f3ffe  0xd11c53907b647e8c  0x7ab22bb43b15585c )
hash ( 0x6f159dfd0835443d  0x55d8a9877b672439  0x327dc934599e4bb7  0xa50a8c9661151fe5  0x3b1224409e600901  0x2f2e2d34681eb87b  0x28de89ca8e928bc2  0xd3ec419ff5954d92  0  0  0  0 ) -> 
        ( 0xfc0d5425584071af  0xffeb54efc8fc83c8  0xf0f2465d67b66762  0xaf7bf91c1565eb2f )
##########
hash ( 0x2e6407ab14535ea  0xa07ff473834f3ffe  0xd11c53907b647e8c  0x7ab22bb43b15585c  0xfc0d5425584071af  0xffeb54efc8fc83c8  0xf0f2465d67b66762  0xaf7bf91c1565eb2f  0  0  0  0 ) -> 
        ( 0x6790b537f607e9b1  0x6a4061557f6e9ae2  0x9008ee29d0ead7b  0xd5d920fb1c6e87ba )
```