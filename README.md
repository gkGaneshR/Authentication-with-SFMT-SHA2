# Authentication-with-SFMT-SHA2
The aim for this  work  is to provide an improved method of  the  S/KEY  while  maintaining  the  balance  between efficiency and protection towards attacks

The password input by the user is given as seed to SIMD-oriented Fast Mersenne Twister (SFMT) to generate
random number. Thus, each bit is aligned with the corresponding first twisted value. If the bit is ‘1’, 
then the first twisted value is again twisted to obtain second twisted value, if ‘0’, then the first twisted 
value is maintained. The single concatenated value is hashed using SHA-2. The hashed value is then input into
SFMT once more and obtain a single value/ password, OTP. User then send the OTPN-user to server, server
compares the obtained password with the one it has, OTPN-server. If the OTPN obtained by server is same as 
the one with user, OTPN-user = OTPNserver, then the password is authenticated and the user gain access. 
Else, the connection is terminated. After 2^20 logins has been completed (after 1, 048, 575 logins has been
done), the system will reinitialised to the initial stage and user will be prompt to change the password.
