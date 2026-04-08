#!/usr/bin/env perl
# A modified rest_auth.pl by Dr. Xiangyang Li in BIGSdb (https://github.com/kjolley/BIGSdb) to download data from PubMLST via REST interface with automatic token renewal 

#The test databases can be reached at https://pubmlst.org/test/.
#To use these, sign up for a PubMLST account (https://pubmlst.org/site_accounts.shtml)
#and link this account with the pubmlst_test_seqdef and pubmlst_test_isolates
#databases (https://pubmlst.org/site_accounts.shtml#registering_with_databases)
#Version 20250130

use strict;
use warnings;
use 5.010;
use FindBin;
use lib "$FindBin::Bin/lib";
use Net::OAuth 0.20;
$Net::OAuth::PROTOCOL_VERSION = Net::OAuth::PROTOCOL_VERSION_1_0A;
use HTTP::Request::Common;
use LWP::UserAgent;
use JSON qw(encode_json decode_json);
use Data::Random qw(rand_chars);
use Data::Dumper;
use Config::Tiny;
use Getopt::Long qw(:config no_ignore_case);
use Term::Cap;
use POSIX;
use File::Basename;

# OAuth credentials for PubMLST
use constant CONSUMER_KEY    => 'A0IPaxPdVT84GyXeeXG8Iup2';
use constant CONSUMER_SECRET => 'fgdmRPXCAHFOaMrghi2wZe57G3JkuIfiTOOcVxP2xL';
use constant REST_BASE_URL   => 'https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef';

my %opts;
GetOptions(
    'a|arguments=s' => \$opts{'a'},    # Additional query parameters
    'm|method=s'    => \$opts{'m'},    # HTTP method (default GET)
    'r|route=s'     => \$opts{'r'},    # API endpoint route
    'o|output=s'    => \$opts{'o'},    # Output file path
    'h|help'        => \$opts{'h'},    # Help message
) or die("Error in command line arguments\n");

if ($opts{'h'}) {
    show_help();
    exit;
}

# Validate HTTP method
$opts{'m'} //= 'GET';
die "Only GET method is supported for downloads.\n" if $opts{'m'} ne 'GET';

main();

sub main {
    my $route = $opts{'r'} || die "Route parameter is required (e.g. --route 'loci')\n";
    
    my ($session_token, $session_secret);
    my $data;
    my $retry_count = 0;
    
    while ($retry_count < 2) {  # Allow one retry
        # Retrieve or create session token
        ($session_token, $session_secret) = _retrieve_token('session_token');
        if (!defined $session_token || !defined $session_secret) {
            my $session_response = _get_session_token();
            ($session_token, $session_secret) = ($session_response->token, $session_response->token_secret);
        }
        
        # Attempt to download data
        my $response = _download_data($route, $session_token, $session_secret);
        
        if ($response->is_success) {
            $data = $response->content;
            last;
        } elsif ($response->code == 401) {  # Unauthorized, token likely expired
            warn "Session token expired or invalid. Renewing...\n";
            unlink 'session_token' if -e 'session_token';
            ($session_token, $session_secret) = (undef, undef);
            $retry_count++;
        } else {
            die "Download failed: " . $response->status_line . "\n";
        }
    }
    
    die "Failed to download data after retry. Check credentials or network.\n" unless defined $data;
    
    # Save or output data
    if ($opts{'o'}) {
        open(my $fh, '>', $opts{'o'}) or die "Cannot open output file: $!";
        print $fh $data;
        close $fh;
        #say "Data saved to: $opts{'o'}";
    } else {
        print $data;
    }
}

sub _download_data {
    my ($route, $session_token, $session_secret) = @_;
    
    my $url = REST_BASE_URL . "/$route";
    
    # Add query parameters
    my %params;
    if ($opts{'a'}) {
        my @pairs = split /&/, $opts{'a'};
        for my $pair (@pairs) {
            my ($key, $value) = split /=/, $pair;
            $params{$key} = $value;
        }
    }
    
    #say "Downloading data from: $url";
    
    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $session_token,
        token_secret     => $session_secret,
        request_url      => $url,
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        extra_params     => \%params,
    );
    
    $request->sign;
    die "OAuth signature verification failed!\n" unless $request->verify;
    
    my $ua = LWP::UserAgent->new;
    return $ua->get($request->to_url);
}

sub _get_request_token {
    my $request = Net::OAuth->request('request token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        request_url      => REST_BASE_URL . '/oauth/get_request_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        callback         => 'oob'
    );
    $request->sign;
    die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
    
    say 'Getting request token...';
    my $ua = LWP::UserAgent->new;
    my $res = $ua->request(GET $request->to_url);
    
    if ($res->is_success) {
        my $decoded_json = decode_json($res->content);
        my $request_response = Net::OAuth->response('request token')->from_hash($decoded_json);
        say 'Request Token:        ', $request_response->token;
        say 'Request Token Secret: ', $request_response->token_secret;
        _write_token('request_token', $request_response->token, $request_response->token_secret);
        return $request_response;
    } else {
        die "Failed to get request token: " . $res->status_line . "\n";
    }
}

sub _get_access_token {
    my ($request_token, $request_secret) = @_;
    
    unless ($request_token && $request_secret) {
        ($request_token, $request_secret) = _retrieve_token('request_token');
        unless ($request_token && $request_secret) {
            my $request_response = _get_request_token();
            ($request_token, $request_secret) = ($request_response->token, $request_response->token_secret);
        }
    }
    
    say "\nPlease visit the following URL to authorize:";
    say "https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=authorizeClient&oauth_token=$request_token";
    print "\nEnter verification code: ";
    my $verifier = <STDIN>;
    chomp $verifier;
    
    my $request = Net::OAuth->request('access token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $request_token,
        token_secret     => $request_secret,
        verifier         => $verifier,
        request_url      => REST_BASE_URL . '/oauth/get_access_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    
    $request->sign;
    die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
    
    say "\nExchanging for access token...";
    unlink 'request_token' if -e 'request_token';
    
    my $ua = LWP::UserAgent->new;
    my $res = $ua->request(GET $request->to_url);
    
    if ($res->is_success) {
        my $decoded_json = decode_json($res->content);
        my $access_response = Net::OAuth->response('access token')->from_hash($decoded_json);
        say 'Access Token:        ', $access_response->token;
        say 'Access Token Secret: ', $access_response->token_secret;
        _write_token('access_token', $access_response->token, $access_response->token_secret);
        return $access_response;
    } else {
        die "Failed to get access token: " . $res->status_line . "\n";
    }
}

sub _get_session_token {
    my ($access_token, $access_secret) = @_;
    
    unless ($access_token && $access_secret) {
        ($access_token, $access_secret) = _retrieve_token('access_token');
        unless ($access_token && $access_secret) {
            my $access_response = _get_access_token();
            ($access_token, $access_secret) = ($access_response->token, $access_response->token_secret);
        }
    }
    
    say "\nRequesting session token...";
    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $access_token,
        token_secret     => $access_secret,
        request_url      => REST_BASE_URL . '/oauth/get_session_token',
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    
    $request->sign;
    die "COULDN'T VERIFY! Check OAuth parameters.\n" unless $request->verify;
    
    my $ua = LWP::UserAgent->new;
    my $res = $ua->request(GET $request->to_url);
    
    if ($res->is_success) {
        my $decoded_json = decode_json($res->content);
        my $session_response = Net::OAuth->response('access token')->from_hash($decoded_json);
        say 'Session Token:        ', $session_response->token;
        say 'Session Token Secret: ', $session_response->token_secret;
        _write_token('session_token', $session_response->token, $session_response->token_secret);
        return $session_response;
    } else {
        die "Failed to get session token: " . $res->status_line . "\n";
    }
}

sub _retrieve_token {
    my ($token_name) = @_;
    return unless -e $token_name;
    
    my $config = Config::Tiny->read($token_name) or 
        die "Failed to read token file: $token_name\n";
    
    return ($config->{_}->{token}, $config->{_}->{secret});
}

sub _write_token {
    my ($token_name, $token, $secret) = @_;
    
    my $config = Config::Tiny->new;
    $config->{_} = {
        token  => $token,
        secret => $secret,
    };
    
    $config->write($token_name) or 
        die "Failed to write token file: $token_name\n";
    
    return 1;
}

sub show_help {
    print << "HELP";
PubMLST Data Download Tool

Usage:
  $0 --route ENDPOINT [options]

Required:
  --route, -r    API endpoint (e.g. 'loci', 'alleles/abcZ', 'sequence/123')

Options:
  --output, -o   Save output to file (default: print to STDOUT)
  --arguments, -a  Query parameters (e.g. 'limit=100&offset=0')
  --help, -h     Show this help message

Examples:
  # Download list of loci
  $0 -r loci -o loci.json
  
  # Download alleles for abcZ locus
  $0 -r alleles/abcZ -o abcZ_alleles.fasta
  
  # Download specific sequence with pagination
  $0 -r sequences -a 'limit=50&offset=100' -o sequences_page2.json

Authentication:
  On first run, you'll be guided through OAuth authentication process.
  Subsequent runs will reuse stored credentials until they expire.

HELP
    exit;
}
