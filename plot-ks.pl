#!/usr/bin/perl
use strict;
use warnings;
use POSIX;
use IO::Pipe;
use IO::Select;
use Getopt::Long;
use Cwd qw(abs_path);
use File::Path qw(remove_tree);
use Storable qw(fd_retrieve store_fd);

# Most processes we run simultaneously
my $max_forks = get_free_cpus();

# Store parent pid
my $parent_pid = $$;

# Hash used to perform reverse translations
my %rev_codon_table = (
	S => qr/((AG[C|T])|(TC.))/,
	F => qr/TT[T|C]/,
	L => qr/((TT[A|G])|(CT.))/,
	Y => qr/TA[T|C]/,
	C => qr/TG[T|C]/,
	W => qr/TGG/,
	P => qr/CC./,
	H => qr/CA[T|C]/,
	Q => qr/CA[A|G]/,
	R => qr/((AG[A|G])|(CG.))/,
	I => qr/AT[T|C|A]/,
	M => qr/ATG/,
	T => qr/AC./,
	N => qr/AA[T|C]/,
	K => qr/AA[A|G]/,
	V => qr/GT./,
	A => qr/GC./,
	D => qr/GA[T|C]/,
	E => qr/GA[A|G]/,
	G => qr/GG./,
	X => qr/.../,
);

# Pfam settings
my $pfam_search;
my $pfam_cpus = $max_forks;

# Allowed models 
my %models = (NG => 1, LWL => 1, LPB => 1, MLWL => 1, MLPB => 1, 
              GY => 1, YN => 1, MYN => 1, MS => 1, MA => 1);

# Default settings
my $model = "YN";
my $bin_size = 0.05;
my $match_length_threshold = 300; # Nucleotides

# Range for Ks plot
my $ks_min = 0;
my $ks_max = 3;
my $exclude_zero;

# Allow for reusing info from an old run
my $input_is_dir = 0;

# Name of output directory
my $project_name = "plot-ks-".time();

my $kaks_calc_in_dir = "kaks-in";
my $kaks_calc_out_dir = "kaks-out";
my $split_blat_out_dir = "blat-out";

# Check that dependencies can be found in user's PATH
my $r = check_path_for_exec("R");
my $blat = check_path_for_exec("blat");
my $kaks_calculator = check_path_for_exec("KaKs_Calculator");

# Check that a version of Transdecoder is present in user's PATH
my $transdecoder = check_path_for_exec("TransDecoder", 1);
my $transdecoder_orfs = check_path_for_exec("TransDecoder.LongOrfs", 1);
my $transdecoder_predict = check_path_for_exec("TransDecoder.Predict", 1);
if (!defined($transdecoder) || (!defined($transdecoder_orfs) && !defined($transdecoder_predict))) {;
	die "Could not locate required TransDecoder executables (TransDecoder or TransDecoder.LongOrfs and TransDecoder.Predict) in PATH.\n" 
}

# Save how script was invoked
my $settings = "@ARGV";

# Parse command line options
GetOptions(
	"model|m=s" => \$model,
	"ks-min=f" => \$ks_min,
	"ks-max=f" => \$ks_max,
	"bin-size|b=f" => \$bin_size,
	"out-dir=s" => \$project_name,
	"exclude-zero|x" => \$exclude_zero,
	"min-length|l=i" => \$match_length_threshold,
	"n-threads|T=i" => \$max_forks,
	"pfam-cpus|p=i" => \$pfam_cpus,
	"pfam-search|s=s" => \$pfam_search,
	"help|h" => \&help,
	"usage" => \&usage
);

# Check that hmmscan can be located if user wants pfam search
my $hmmscan = check_path_for_exec("hmmscan") if ($pfam_search);

# Error check input
my $transcriptome = shift;
die "You must specify a FASTA file containing transcriptome data for input OR the output directory from the previous execution of this script.\n".&usage if (!defined($transcriptome));
die "Could not locate '$transcriptome'.\n".&usage if (!-e $transcriptome);
die "Your Ks plot bin size can't be a negative number.\n".&usage if ($bin_size < 0);
die "Your Ks plot bin size is larger than your Ks range.\n".&usage if ($bin_size > $ks_max - $ks_min);
die "The model '$model' does not exist or is not usable by this script.\n".&usage if (!exists($models{$model}));
die "Could not locate hmmscan binary file: '$pfam_search'.\n" if (defined($pfam_search) && !-e $pfam_search);

# Make sure the TransDecoder pfam script can be found if needed
check_path_for_exec("pfam_runner.pl") if (defined($pfam_search));

# Check if the ks span is evenly divisible by the bin size
if ($bin_size != 0 && (($ks_max - $ks_min) / $bin_size) =~ /(\.\d+)/) {
	my $remainder = $1 * $bin_size;
	logger("WARNING: The span of Ks values you wished to plot ($ks_min to $ks_max) was not evenly divisible by your bin size ($bin_size).");
	$ks_max -= $remainder;
	logger("The maximum possible Ks value which will be plotted has therefore been reduced to $ks_max.");
}

undef($exclude_zero) if ($ks_min > 0);

# Input is a previous run directory, reuse information
$input_is_dir++ if (-d $transcriptome);

my $input_root;
if (!$input_is_dir) {

	# Clean run with no prior output
	my $transcriptome_abs_path = abs_path($transcriptome);

	# Initialize our working directory and create a symlink to input
	mkdir($project_name) if (!-e $project_name) || die "Could not create directory '$project_name': $!.\n";

	$transcriptome =~ s/.*\///;
	run_cmd("ln -s $transcriptome_abs_path $ENV{PWD}/$project_name/$transcriptome");
	($input_root = $transcriptome) =~ s/(.*)\.\S+$/$1/;
}
else {

	# Prior output available, set relevant variables

	$project_name = $transcriptome;
	my @contents = glob("$project_name/*");

	# Determine the transcriptome name by looking for a symlink
	my $found_name = 0;
	foreach my $file (@contents) {
		if (-l $file) {
			$file =~ s/\Q$project_name\E\///;
			$transcriptome = $file;
			$found_name = 1;
		}
	}
	($input_root = $transcriptome) =~ s/(.*)\.\S+$/$1/;

	die "Could not locate transcriptome in '$project_name'.\n" if (!$found_name);
}

chdir($project_name);

# Change how we handle Ctrl+C
$SIG{INT} = 'INT_handler';

# Set up output file names
my $blat_out = "self-blat.pslx";
my $paralogs_seqs = "paralogs-t$match_length_threshold.atx";
my $paralogs_ka_ks_out = "paralogs-t$match_length_threshold-m$model.kaks";
my $final_ks_values = "paralogs-t$match_length_threshold-m$model.csv";

# Output plot name varies based on bin size and Ks range
my $ks_plot_name = "$input_root-t$match_length_threshold-m$model";
if ($exclude_zero) {
	if ($bin_size == 0) {
		$ks_plot_name .= "-density-range\($ks_min-$ks_max].pdf";
	}
	else {
		$ks_plot_name .= "-b$bin_size-range\($ks_min-$ks_max].pdf";
	}
}
else {
	if ($bin_size == 0) {
		$ks_plot_name .= "-density-range[$ks_min-$ks_max].pdf";
	}
	else {
		$ks_plot_name .= "-b$bin_size-range[$ks_min-$ks_max].pdf";
	}
}

my $pid = fork();

# Child fork
if ($pid == 0) {

	# Cleans up if cancelled
	$SIG{TERM} = 'TERM_handler';

	# Echo script invocation
	logger("Invocation: perl plot-ks.pl $settings");

	# Translate transcriptome into most likely proteome
	run_transdecoder();

	# Identify paralogs
	run_self_blat();
	parse_self_blat_output();

	# Calculate Ks for paralogs
	run_kaks_calc();
	parse_ks_values();

	# Create the final plot
	create_ks_plot();
	exit(0);
}
else {
	waitpid($pid, 0);
}

sub run_transdecoder {

	my $transdecoder_out_dir = "transdecoder";

	# Check if the files we need from TransDecoder are present
	if (!-e "$transcriptome.transdecoder.pep" || !-e "$transcriptome.transdecoder.mRNA") {
		logger("\nRunning TransDecoder on '$transcriptome'...");

		# Opt for newer version of TransDecoder if present
		if ($transdecoder_orfs) {
			run_cmd("$transdecoder_orfs -t $transcriptome");

			# Need to run things differently to include pfam
			if (defined($pfam_search)) {
				run_cmd("$hmmscan --cpu $pfam_cpus --domtblout pfam.domtblout $pfam_search $transcriptome.transdecoder_dir/longest_orfs.pep");
				run_cmd("$transdecoder_predict -t $transcriptome --retain_pfam_hits pfam.domtblout");
			}
			else {
				run_cmd("$transdecoder_predict -t $transcriptome");
			}
			remove_tree("$transcriptome.transdecoder_dir/");
		}
		# Old version of TransDecoder
		else {
			if (defined($pfam_search)) {
				run_cmd("$transdecoder -t $transcriptome --workdir $transdecoder_out_dir --search_pfam $pfam_search --CPU $pfam_cpus");
			}
			else {
				run_cmd("$transdecoder -t $transcriptome --workdir $transdecoder_out_dir");
			}
			remove_tree($transdecoder_out_dir);
		}

		# Clean up unneeded files
		remove_tree($transdecoder_out_dir);

		logger("Completed TransDecoder.\n");
	}
	else {
		logger("\nUsing TransDecoder output from previous run.");
	}
}

sub run_self_blat {

	# Check if we've already run self-blat
	if (!-e $blat_out) {
		logger("Running self-blat...");
		run_cmd("$blat $transcriptome.transdecoder.pep $transcriptome.transdecoder.pep -prot -out=pslx self-blat.pslx -noHead");
		logger("Completed self-blat.\n");
	}
	else {
		logger("Using blat output from previous run located in '$blat_out'.");
	}
}

sub parse_self_blat_output {

	# Check if we've already parsed self-blat output
	if (!-e $paralogs_seqs) {

		logger("Parsing self-blat output...");

		# Determine how to split the input file
		chomp(my $total_hits = `wc -l $blat_out`);
		$total_hits =~ s/(^\s*\d+).*/$1/;

		my $lines_per_thread = ceil($total_hits / $max_forks);

		# Create the directory that will hold the split blat output
		mkdir("$split_blat_out_dir") || die "Could not create directory '$split_blat_out_dir': $!.\n";

		# Split the output from blat
		run_cmd("split -l $lines_per_thread $blat_out $split_blat_out_dir/$blat_out.");

		# Run each parition of the blat output concurrently

		my @pids;
		my $select = new IO::Select;
		foreach my $file (glob("$split_blat_out_dir/$blat_out.*")) {

			# We need a pipe to send hashes from children to parent
			my $pipe = new IO::Pipe;	
			my $pid = fork();

			if ($pid == 0) {
				# Child
				my $to_parent = $pipe->writer();
				parse_self_blat_output_child($file, $to_parent);

				# Tell the Parent we have finished
				until(Storable::is_storing() == 0){};
				store_fd({DONE => $$}, $to_parent);

				exit(0);
			}
			else {
				# Parent
				my $from_child = $pipe->reader();
				$select->add($from_child);
				push(@pids, $pid);
			}
		}

		# Retrieve and parse output returned by forks

		my %queries;
		my %matches;
		while (my @handles = $select->can_read) {

			my $handle = shift(@handles);
			my %hit = %{fd_retrieve($handle)};

			# Fork has parsed output, close and remove pipe to it
			if (exists($hit{"DONE"})) {
				waitpid($hit{"DONE"}, 0);
				$select->remove($handle);
				$handle->close();
				next;
			}
			else {
				my $query_name = $hit{'QUERY_NAME'};
				my $match_name = $hit{'MATCH_NAME'};

				# Check for hit with same pair but with query/match reversed
				if (!exists($queries{"$match_name-$query_name"})) {

					my $name = "q_$query_name"."_t_$match_name";

					# Check if there is already a match between these two seque
					# if there is a match, the longer length one will be output
					if (exists($matches{$name})) {
						my $current_length = $matches{$name}->{'LENGTH'};
						if ($current_length <= $hit{'LENGTH'}) {
							$matches{$name} = \%hit;
						}
					}
					else {
						$matches{$name} = \%hit;
					}
					
					# Keep track of pairs so we can tell if a pair is
					# the same but with query/match reversed
					$queries{"$query_name-$match_name"}++;
				}
			}
		}

		# Just for good measure
		foreach my $pid (@pids) {
			waitpid($pid, 0);
		}

		# Total number of pairs found
		my $final_hits = scalar(keys %matches);
		
		# Ensure that we actually got hits
		if ($final_hits == 0) {
			logger("No blat hits met the requirements.\n");
			exit(0);
		}
		logger("Completed parsing blat output, $final_hits blat hit(s) met the requirements.\n");

		# Output parsed hits	
		open(my $output_file, ">", $paralogs_seqs);
		foreach my $key (sort { $a cmp $b} keys %matches) {
			my $pair = $matches{$key};

			print {$output_file} ">$key\n";
			print {$output_file} "$pair->{QUERY_ALIGN}\n";
			print {$output_file} "$pair->{MATCH_ALIGN}\n\n";
		}
		close($output_file);

		# Clean up split files
		remove_tree($split_blat_out_dir);
	}
	else {
		logger("Using previously parsed blat output in '$paralogs_seqs'.");
	}
}

sub parse_self_blat_output_child {
	my ($file, $handle) = @_;

	# Load transcriptome sequences
	my %align = parse_fasta("$transcriptome.transdecoder.mRNA");

	# Parse this fork's portion of the blat output
	my %matches;
	open(my $blat_out, "<", $file);
	while (my $line = <$blat_out>) {
		chomp($line);

		# Split the line on tabs and extract the information we want
		my @line = split("\t", $line);
		my ($query_name, $match_name, $query_align, $match_align) = 
			($line[9], $line[13], $line[21], $line[22]);

		# Don't want self matches
		next if ($query_name eq $match_name);
		
		# Reverse translate amino acids to their corresponding nucleotides

		my @query_align =  split(",", $query_align);
		my @match_align =  split(",", $match_align);

		my $query_nuc_align = $align{$query_name};
		my $match_nuc_align = $align{$match_name};

		my $trans_query_align;
		foreach my $align (@query_align) {
			$trans_query_align .= reverse_translate({"DNA" => $query_nuc_align, "PROT" => $align});
		}

		my $trans_match_align;
		foreach my $align (@match_align) {
			$trans_match_align .= reverse_translate({"DNA" => $match_nuc_align, "PROT" => $align});
		}

		$query_align = $trans_query_align;
		$match_align = $trans_match_align;

		# Check that length threshold is met
		if (length($query_align) >= $match_length_threshold) {
			
			# Create a hash containing all needed information on this hit
			my $pair = {'QUERY_NAME'  => $query_name,
						'MATCH_NAME'  => $match_name,
						'QUERY_ALIGN' => $query_align,
						'MATCH_ALIGN' => $match_align,
						'LENGTH' => length($query_align)};

			# Send the hash back to the parent process
			until(Storable::is_storing() == 0){};
			store_fd($pair, $handle);
		}
	}
}


sub run_kaks_calc {

	# Check if we've run KaKs_Calculator before
	if (!-e $paralogs_ka_ks_out) {

		# Run KaKs_Calculator to get Ks values
		logger("Running KaKs_Calculator...");

		# Determine how to split the input file
		chomp(my $total_hits = `wc -l $paralogs_seqs`);
		$total_hits =~ s/(^\s*\d+).*/$1/;
		$total_hits = $total_hits / 4;

		my $lines_per_thread = ceil($total_hits / $max_forks) * 4;

		# Create the directories that will hold the split files and their KaKs_Calculator output
		mkdir("$kaks_calc_in_dir") || die "Could not create directory '$kaks_calc_in_dir': $!.\n";
		mkdir("$kaks_calc_out_dir") || die "Could not create directory '$kaks_calc_out_dir': $!.\n";

		# Split the files
		run_cmd("split -l $lines_per_thread $paralogs_seqs $kaks_calc_in_dir/$paralogs_seqs.");

		# Run each instance of KaKs_Calculator

		my @pids;
		foreach my $file (glob("$kaks_calc_in_dir/$paralogs_seqs.*")) {

			(my $file_id = $file) =~ s/.*\.(\S+$)/$1/;

			my $pid = fork();
			if ($pid == 0) {
				# Child 
				run_cmd("$kaks_calculator -i $file -o $kaks_calc_out_dir/$paralogs_ka_ks_out.$file_id -m $model >/dev/null");
				exit(0);
			}
			else {
				# Parent
				push(@pids, $pid);
			}
		}

		# Don't let forks become zombies
		foreach my $pid (@pids) {
			waitpid($pid, 0);
		}

		# Combine all the output files generated by the multiple forks
		run_cmd("cat $kaks_calc_out_dir/$paralogs_ka_ks_out.* > $paralogs_ka_ks_out");

		# Clean up
		remove_tree($kaks_calc_in_dir);
		remove_tree($kaks_calc_out_dir);

		logger("Completed KaKs_Calculator.\n");
	}
	else {
		logger("Using KaKs_Calculator output from previous run in '$paralogs_ka_ks_out'.");
	}
}

sub parse_ks_values {

	# Check if we've parsed Ks values before
	if (!-e $final_ks_values) {

		# Open KaKs_Calculator output, parse Ks values and convert to .csv for R
		open(my $ka_ks_calculator_output, "<", $paralogs_ka_ks_out);
		open(my $ks_csv, ">", $final_ks_values);
		while (my $line = <$ka_ks_calculator_output>) {
			chomp($line);
			my @line = split(/\s+/, $line);

			# Skip first line containing headers
			if ($. == 1) {
				next;
			}

			next if ($line[0] eq "Sequence");

			my $ks = $line[3];
			$ks = 0 if ($ks eq "NA");

			print {$ks_csv} "$ks\n"
		}
		close($ks_csv);
		close($ka_ks_calculator_output);
	}
}

sub create_ks_plot {

	# Create PDF plot of output
	logger("Creating Ks plot in R...");
	if ($exclude_zero) {
		if ($bin_size == 0) {
			run_cmd("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values', header=F); 
				data <- data\\\$V1[data\\\$V1 <= $ks_max & data\\\$V1 > 0]; 
				plot(density(data), main=expression(paste('K'[s], ' Density Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
		else {
			run_cmd("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values', header=F); 
				data <- data\\\$V1[data\\\$V1 <= $ks_max & data\\\$V1 > 0]; 
				hist(data, breaks=seq($ks_min,$ks_max,by=$bin_size), 
					main=expression(paste('K'[s], ' Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
	}
	else {
		if ($bin_size == 0) {
			run_cmd("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values', header=F); 
				data <- data\\\$V1[data\\\$V1 <= $ks_max]; 
				data <- data[data >= $ks_min];
				plot(density(data), main=expression(paste('K'[s], ' Density Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
		else {
			run_cmd("echo \"pdf(file='$ks_plot_name'); 
				data=read.csv('$final_ks_values', header=F); 
				data <- data\\\$V1[data\\\$V1 <= $ks_max]; 
				data <- data[data >= $ks_min];
				hist(data, breaks=seq($ks_min,$ks_max,by=$bin_size), 
					main=expression(paste('K'[s], ' Plot for $transcriptome')), 
					xlab=expression(paste('Pairwise', ' K'[s])), axes=T);\" | $r --no-save") && die;
		}
	}
	logger("\nCompleted Ks plot.\n");
}

sub reverse_translate {
	my $settings = shift;

	my $dna = $settings->{'DNA'};
	my $prot = $settings->{'PROT'};

	# Generate Regex to perform reverse translation

	my $regex;
	foreach my $index (0 .. length($prot) - 1) {
		my $char = substr($prot, $index, 1);
		$regex .= $rev_codon_table{$char};
	}

	# Find the match of the Regex

	my $translation;
	if ($dna =~ /($regex)/) {
		$translation = $1;
	}
	elsif (reverse($dna) =~ /($regex)/) {
		$translation = $1;
	}
	else {
		die "Protein sequence could not be reverse translated.\n";
	}

	return $translation;
}

sub check_path_for_exec {
	my ($exec, $continue) = @_;
	
	my $path = $ENV{PATH}.":."; # include current directory as well
	my @path_dirs = split(":", $path);

	my $exec_path;
	foreach my $dir (@path_dirs) {
		$dir .= "/" if ($dir !~ /\/$/);
		$exec_path = $dir.$exec if (-e $dir.$exec && -x $dir.$exec && !-d $dir.$exec);
	}

	if (!defined($exec_path) && $continue) {
		return;
	}
	elsif (defined($exec_path)) {
		return $exec_path;
	}
	else {
		die "Could not find the following executable: '$exec'. This script requires this program in your path.\n";
	}
	return $exec_path;
}

sub parse_fasta {
	my $filename = shift;

	my $taxon;
	my %align;
	open(my $alignment_file, '<', $filename) 
		or die "Could not open '$filename': $!\n";

	while (my $line = <$alignment_file>) {
		$line =~ s/^\s+|\s+$//g;

		# Taxon name
		if ($line =~ /^>(\S+)/) {
			$taxon = $1;
		}
		else {
			# Taxon sequence
			$taxon =~ s/-/_/g;
			$align{$taxon} .= $line;
		}
	}
	close($alignment_file);
	
	return %align;
}

sub logger {
	my $msg = shift;

	my $time = "[".localtime(time())."]";

	# Allow for new lines and carriage returns before message
	if ($msg =~ s/^\n//) {
		$time = "\n$time";
	}
	elsif ($msg =~ s/^\r//) {
		$time = "\r$time";
	}

	print "$time $msg\n"; 
}

sub INT_handler {
	kill(15, $pid);
	exit(0);
}

sub TERM_handler {

	logger("\rKeyboard interrupt detected, stopping analyses and cleaning up.");

	# We only want to clean up if the initial input was a sequence file
	if (!$input_is_dir) {

		# Move back into the directory script was called in
		chdir("..");

		# Try to delete directory five times, if it can't be deleted print an error message
		my $count = 0;
		until (!-e $project_name || $count == 5) {
			$count++;

			remove_tree($project_name, {error => \my $err});
			sleep(1);
		}
		logger("Could not clean all files in './$project_name/'.") if ($count == 5);
	}
	else {
		# Clean up multithreading directories
		my $count = 0;
		until (!-e $kaks_calc_in_dir || $count == 5) {
			$count++;

			remove_tree($kaks_calc_in_dir, {error => \my $err});
			sleep(1);
		}
		logger("Could not clean all files in './$project_name/'.") if ($count == 5);
		undef($count);

		until (!-e $kaks_calc_out_dir || $count == 5) {
			$count++;

			remove_tree($kaks_calc_out_dir, {error => \my $err});
			sleep(1);
		}
		logger("Could not clean all files in './$project_name/'.") if ($count == 5);
		undef($count);

		until (!-e $split_blat_out_dir || $count == 5) {
			$count++;

			remove_tree($split_blat_out_dir, {error => \my $err});
			sleep(1);
		}
		logger("Could not clean all files in './$project_name/'.") if ($count == 5);
		undef($count);
	}

	exit(0);
}

sub get_free_cpus {

	my $os_name = $^O;

	# Returns a two-member array containing CPU usage observed by top,
	# top is run twice as its first output is usually inaccurate
	my @percent_free_cpu;
	if ($os_name eq "darwin") {
		# Mac OS
		chomp(@percent_free_cpu = `top -i 1 -l 2 | grep "CPU usage"`);
	}
	else {
		# Linux
		chomp(@percent_free_cpu = `top -b -n2 -d0.05 | grep "Cpu(s)"`);
	}

	my $percent_free_cpu = pop(@percent_free_cpu);

	if ($os_name eq "darwin") {
		# Mac OS
		$percent_free_cpu =~ s/.*?(\d+\.\d+)%\s+id.*/$1/;
	}
	else {
		# linux 
		$percent_free_cpu =~ s/.*?(\d+\.\d)\s*%?ni,\s*(\d+\.\d)\s*%?id.*/$1 + $2/; # also includes %nice as free 
		$percent_free_cpu = eval($percent_free_cpu);
	}

	my $total_cpus;
	if ($os_name eq "darwin") {
		# Mac OS
		$total_cpus = `sysctl -n hw.ncpu`;
	}
	else {
		# linux
		$total_cpus = `grep --count 'cpu' /proc/stat` - 1;
	}

	my $free_cpus = ceil($total_cpus * $percent_free_cpu / 100);

	if ($free_cpus == 0 || $free_cpus !~ /^\d+$/) {
		$free_cpus = 1; # assume that at least one cpu can be used
	}
	
	return $free_cpus;
}

sub run_cmd {
	my $command = shift;

	my $return = system($command);

	if ($return) {
		logger("'$command' died with error: '$return'.\n");
		kill(2, $parent_pid);
		exit(0);
	}
}

sub usage {
	return "Usage: perl plot-ks.pl [TRANSCRIPTOME] [OPTIONS]...\n";
}

sub help {
print <<EOF; 
@{[usage()]}
Generate a Ks plot for a given transcriptome in fasta format

  -m, --model                       model used by KaKs_Calculator to determine Ks (default: YN)
                                        valid models: NG, LWL, LPB, MLWL, MLPB, GY, YN, MYN, MS, MA. See KaKs_Calculator doc for details.
  -l, --min-length                  the minimum alignment length of paralogous sequences (default: 300 bp)
  -x, --exclude-zero                used to exclude Ks = 0 from plot, useful for unfiltered Trinity assemblies
  -b, --bin-size                    size of bins used in histogram of Ks plot, set to 0 for a density plot (default: 0.05)
  -o, --out-dir                     name of the directory to store output files in (default: "plot-ks" + Unix time of script invocation)
  --ks-min                          lower boundary for x-axis of Ks plot (default: Ks = 0)
  --ks-max                          upper boundary for x-axis of Ks plot (default: Ks = 3)
  --pfam-cpus                       the number of CPUs to let hmmscan using during pfam analysis (default: current number of free CPUs)
  --pfam-search                     full path to pfam binary for usage in pfam search
  -T, --n-threads                   the number of CPUs to use during analysis (default: current number of free CPUs)
  -h, --help                        display this help and exit

Examples:
  perl plot-ks.pl assembly.fa -b 0.01               generates a Ks (YN model) plot from [0, 3] using a bin size of 0.01, and contigs 
                                                    with at least 300 bp of homologous sequence

  perl plot-ks.pl assembly.fa -x -m NG              generates a Ks (NG model) plot from (0, 3] using a bin size of 0.05, and contigs 
                                                    with at least 300 bp of homologous sequence

  perl plot-ks.pl assembly.fa --ks-max 5 -t 500     generates a Ks (YN model) plot from [0, 5] using a bin size of 0.01, and contigs 
                                                    with at least 500 bp of homologous sequence

Mail bug reports and suggestions to <noah.stenz.github\@gmail.com>
EOF
exit(0);
}
