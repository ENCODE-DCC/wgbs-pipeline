# Defines image for running linting (R installation takes ages due to compiling libraries)
# You should never need to use this directly.
# Builds with docker build -t quay.io/encode-dcc/wgbs-pipeline:circleci-lint - < Dockerfile
FROM circleci/python:3.7.3-stretch

# Installs rust with cargo clippy and rustfmt
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y

# R lints
RUN sudo apt-get install -y r-base && \
    sudo Rscript -e 'install.packages("lintr", repos = "http://cran.us.r-project.org")' && \
    sudo Rscript -e 'install.packages("styler", repos = "http://cran.us.r-project.org")' && \
    sudo Rscript -e 'install.packages("docopt", repos = "http://cran.us.r-project.org")' && \
    sudo Rscript -e 'install.packages("git2r", repos = "http://cran.us.r-project.org")'

# Install shellcheck
RUN curl -s https://storage.googleapis.com/shellcheck/shellcheck-stable.linux.x86_64.tar.xz | sudo tar xJ && sudo mv shellcheck-stable/shellcheck /usr/local/bin
