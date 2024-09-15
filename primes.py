import concurrent.futures
import math

PRIMES = [
    # 112272535095293,
    # 112582705942171,
    # 112272535095293,
    # 115280095190773,
    # 115797848077099,
    # 1099726899285419,
    768614336404564651,
    108086391056891903,
    10657331232548839,
    790738119649411319,
    489133282872437279,
    2305843009213693951,
]

def is_prime(n):
    if n < 2:
        return False
    if n == 2:
        return True
    if n % 2 == 0:
        return False

    sqrt_n = int(math.floor(math.sqrt(n)))
    for i in range(3, sqrt_n + 1, 2):
        if n % i == 0:
            return False
    return True

def main():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for number, prime in zip(PRIMES, executor.map(is_prime, PRIMES)):
            print('%d is prime: %s' % (number, prime))
    # for n in PRIMES:
    #     print('%d is prime: %s' % (n, is_prime(n)))

if __name__ == '__main__':
    main()
