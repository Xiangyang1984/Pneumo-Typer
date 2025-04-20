#!/usr/bin/env perl
# PubMLST数据下载工具（最终稳定版）
use strict;
use warnings;
use 5.010;
use Net::OAuth 0.20;
$Net::OAuth::PROTOCOL_VERSION = Net::OAuth::PROTOCOL_VERSION_1_0A;
use JSON qw(decode_json);
use Data::Random qw(rand_chars);
use Config::Tiny;
use Getopt::Long qw(:config no_ignore_case);
use File::Temp qw(tempfile);
use URI::Escape;

# OAuth凭证配置
use constant CONSUMER_KEY    => 'A0IPaxPdVT84GyXeeXG8Iup2';
use constant CONSUMER_SECRET => 'fgdmRPXCAHFOaMrghi2wZe57G3JkuIfiTOOcVxP2xL';
use constant REST_BASE_URL   => 'https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef';

my %opts;
GetOptions(
    'a|arguments=s' => \$opts{a},    # 附加查询参数
    'm|method=s'    => \$opts{m},    # HTTP方法（默认GET）
    'r|route=s'     => \$opts{r},    # API端点路径
    'o|output=s'    => \$opts{o},    # 输出文件路径
    'h|help'        => \$opts{h},    # 显示帮助
    'v|verbose'     => \$opts{v},    # 详细输出
) or die("命令行参数错误\n");

if ($opts{h}) {
    show_help();
    exit;
}

$opts{m} //= 'GET';
die "仅支持GET方法下载\n" if $opts{m} ne 'GET';

main();

sub main {
    my $route = $opts{r} || die "需要指定API路径（例如--route 'loci'）\n";
    
    my ($session_token, $session_secret) = _get_valid_session_token();
    
    my $data = _download_data($route, $session_token, $session_secret);
    
    if ($opts{o}) {
        open(my $fh, '>', $opts{o}) or die "无法打开输出文件: $!";
        print $fh $data;
        close $fh;
        warn "数据已保存至: $opts{o}\n" if $opts{v};
    } else {
        print $data;
    }
}

sub _download_data {
    my ($route, $token, $secret) = @_;
    
    my $url = REST_BASE_URL . "/$route";
    my %params;
    if ($opts{a}) {
        %params = map { split(/=/, $_, 2) } split(/&/, $opts{a});
    }
    
    warn "正在下载: $url\n" if $opts{v};
    
    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $token,
        token_secret     => $secret,
        request_url      => $url,
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        extra_params     => \%params,
    );
    
    $request->sign;
    die "OAuth签名验证失败!\n" unless $request->verify;
    
    # 构建OAuth头
    my @oauth_params = (
        'oauth_consumer_key="'    . uri_escape(CONSUMER_KEY) . '"',
        'oauth_nonce="'           . uri_escape($request->nonce) . '"',
        'oauth_signature="'       . uri_escape($request->signature) . '"',
        'oauth_signature_method="HMAC-SHA1"',
        'oauth_timestamp="'      . uri_escape($request->timestamp) . '"',
        'oauth_version="1.0"',
    );
    push @oauth_params, 'oauth_token="' . uri_escape($token) . '"' if $token;
    
    my $auth_header = 'Authorization: OAuth ' . join(', ', @oauth_params);
    
    my ($out_fh, $out_file) = tempfile(UNLINK => 1);
    close $out_fh;
    
    my $wget_cmd = "wget --quiet --header '$auth_header' " .
                   "--output-document=" . shell_quote($out_file) . " " .
                   "--timeout=30 " . shell_quote($url);
    
    warn "执行命令: $wget_cmd\n" if $opts{v};
    
    my $exit_status = system($wget_cmd);
    die "下载失败 (退出状态: $exit_status)\n" if $exit_status != 0;
    
    open(my $fh, '<', $out_file) or die "无法读取临时文件: $!";
    my $content = do { local $/; <$fh> };
    close $fh;
    
    if ($content =~ /"error":"invalid_token"/) {
        die "会话令牌无效或已过期\n";
    }
    
    return $content;
}

sub _get_valid_session_token {
    my $max_retries = 3;
    for my $attempt (1..$max_retries) {
        my ($token, $secret) = _retrieve_token('session_token');
        
        if ($token && $secret && _validate_token($token, $secret)) {
            return ($token, $secret);
        }
        
        warn "获取新会话令牌（尝试: $attempt）...\n" if $opts{v};
        my $session = _get_session_token();
        _write_token('session_token', $session->token, $session->token_secret);
    }
    die "无法获取有效会话令牌\n";
}

sub _validate_token {
    my ($token, $secret) = @_;
    
    my $test_url = REST_BASE_URL . '/loci?limit=1';
    my $request = Net::OAuth->request('protected resource')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $token,
        token_secret     => $secret,
        request_url      => $test_url,
        request_method   => 'GET',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    
    $request->sign;
    
    my @oauth_params = (
        'oauth_consumer_key="'    . uri_escape(CONSUMER_KEY) . '"',
        'oauth_nonce="'           . uri_escape($request->nonce) . '"',
        'oauth_signature="'       . uri_escape($request->signature) . '"',
        'oauth_signature_method="HMAC-SHA1"',
        'oauth_timestamp="'      . uri_escape($request->timestamp) . '"',
        'oauth_version="1.0"',
        'oauth_token="'          . uri_escape($token) . '"'
    );
    
    my $auth_header = 'Authorization: OAuth ' . join(', ', @oauth_params);
    
    my $cmd = "wget --spider --header " . shell_quote($auth_header) . " " .
              "--timeout=15 " . shell_quote($test_url) . " >/dev/null 2>&1";
    my $exit_status = system($cmd);
    return $exit_status == 0;
}

sub _get_session_token {
    my ($access_token, $access_secret) = _retrieve_token('access_token');
    
    unless ($access_token) {
        my $access_res = _get_access_token();
        ($access_token, $access_secret) = ($access_res->token, $access_res->token_secret);
    }
    
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
    
    my @oauth_params = (
        'oauth_consumer_key="'    . uri_escape(CONSUMER_KEY) . '"',
        'oauth_nonce="'           . uri_escape($request->nonce) . '"',
        'oauth_signature="'       . uri_escape($request->signature) . '"',
        'oauth_signature_method="HMAC-SHA1"',
        'oauth_timestamp="'      . uri_escape($request->timestamp) . '"',
        'oauth_version="1.0"',
        'oauth_token="'          . uri_escape($access_token) . '"'
    );
    
    my $auth_header = 'Authorization: OAuth ' . join(', ', @oauth_params);
    
    my ($out_fh, $out_file) = tempfile(UNLINK => 1);
    close $out_fh;
    
    my $cmd = "wget --quiet --header " . shell_quote($auth_header) . " " .
              "--output-document=" . shell_quote($out_file) . " " .
              shell_quote($request->to_url);
    
    system($cmd) == 0 or die "获取会话令牌失败\n";
    
    open(my $fh, '<', $out_file) or die $!;
    my $json = decode_json(do { local $/; <$fh> });
    return Net::OAuth->response('access token')->from_hash($json);
}

sub _get_access_token {
    my $request_token = _get_request_token();
    
    print "\n请访问以下URL授权：\n";
    print "https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=authorizeClient&oauth_token=".$request_token->token."\n";
    print "输入验证码：";
    my $verifier = <STDIN>;
    chomp $verifier;
    
    my $request = Net::OAuth->request('access token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        token            => $request_token->token,
        token_secret     => $request_token->token_secret,
        verifier         => $verifier,
        request_url      => REST_BASE_URL . '/oauth/get_access_token',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
    );
    
    $request->sign;
    
    my ($out_fh, $out_file) = tempfile(UNLINK => 1);
    close $out_fh;
    
    my $cmd = "wget --quiet --output-document=" . shell_quote($out_file) . " " .
              shell_quote($request->to_url);
    system($cmd) == 0 or die "获取访问令牌失败\n";
    
    open(my $fh, '<', $out_file) or die $!;
    my $json = decode_json(do { local $/; <$fh> });
    _write_token('access_token', $json->{oauth_token}, $json->{oauth_token_secret});
    return Net::OAuth->response('access token')->from_hash($json);
}

sub _get_request_token {
    my $request = Net::OAuth->request('request token')->new(
        consumer_key     => CONSUMER_KEY,
        consumer_secret  => CONSUMER_SECRET,
        request_url      => REST_BASE_URL . '/oauth/get_request_token',
        signature_method => 'HMAC-SHA1',
        timestamp        => time,
        nonce            => join('', rand_chars(size => 16, set => 'alphanumeric')),
        callback         => 'oob'
    );
    
    $request->sign;
    
    my ($out_fh, $out_file) = tempfile(UNLINK => 1);
    close $out_fh;
    
    my $cmd = "wget --quiet --output-document=" . shell_quote($out_file) . " " .
              shell_quote($request->to_url);
    system($cmd) == 0 or die "获取请求令牌失败\n";
    
    open(my $fh, '<', $out_file) or die $!;
    my $content = do { local $/; <$fh> };
    my $response = Net::OAuth->response('request token')->from_post_body($content);
    _write_token('request_token', $response->token, $response->token_secret);
    return $response;
}

sub _retrieve_token {
    my ($name) = @_;
    return unless -e $name;
    
    my $config = Config::Tiny->read($name) or return;
    return ($config->{_}->{token}, $config->{_}->{secret});
}

sub _write_token {
    my ($name, $token, $secret) = @_;
    my $config = Config::Tiny->new;
    $config->{_} = { token => $token, secret => $secret };
    $config->write($name) or die "无法写入令牌文件: $name";
}

sub shell_quote {
    my ($str) = @_;
    $str =~ s/'/'\\''/g;
    return "'$str'";
}

sub show_help {
    print <<'HELP';
PubMLST数据下载工具（OAuth 1.0a认证）

使用方法:
  perl download_pubmlst.pl --route API路径 [选项]

选项:
  -r, --route     必填 API路径（如 loci/alleles_fasta）
  -o, --output    输出文件路径（默认输出到屏幕）
  -a, --arguments 查询参数（如 'limit=100&offset=0'）
  -v, --verbose   显示详细日志
  -h, --help      显示帮助信息

示例:
  # 下载所有位点列表
  perl download_pubmlst.pl -r loci -o loci.json

  # 下载指定基因的等位基因（FASTA格式）
  perl download_pubmlst.pl -r loci/RecA/alleles_fasta -o RecA.fas

  # 带参数查询
  perl download_pubmlst.pl -r profiles -a 'limit=50&fields=ST,clonal_complex' -o cc.csv

HELP
}
