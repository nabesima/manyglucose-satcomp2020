FROM ubuntu:16.04

#RUN apt-get update && apt-get install -y openssh-server git
RUN apt-get update && apt-get install -y git
#RUN mkdir /var/run/sshd
RUN echo 'root:THEPASSWORDYOUCREATED' | chpasswd
#RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
#RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd
#RUN git clone https://github.com/biotomas/hordesat

#ENV NOTVISIBLE "in users profile"
#RUN echo "export VISIBLE=now" >> /etc/profile


# RUN echo "${USER} ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers
# ENV SSHDIR /root/.ssh
# RUN mkdir -p ${SSHDIR}
# RUN touch ${SSHDIR}/sshd_config
# RUN ssh-keygen -t rsa -f ${SSHDIR}/ssh_host_rsa_key -N ''
# RUN cp ${SSHDIR}/ssh_host_rsa_key.pub ${SSHDIR}/authorized_keys
# RUN cp ${SSHDIR}/ssh_host_rsa_key ${SSHDIR}/id_rsa
# RUN echo " IdentityFile ${SSHDIR}/id_rsa" >> /etc/ssh/ssh_config
# RUN echo "Host *" >> /etc/ssh/ssh_config && echo " StrictHostKeyChecking no" >> /etc/ssh/ssh_config
# RUN chmod -R 600 ${SSHDIR}/* && \
# chown -R ${USER}:${USER} ${SSHDIR}/
# # check if ssh agent is running or not, if not, run
# RUN eval `ssh-agent -s` && ssh-add ${SSHDIR}/id_rsa

RUN apt-get update
RUN apt-get install wget -y
RUN apt-get install unzip
RUN apt-get install build-essential -y
RUN apt-get install zlib1g-dev -y
RUN DEBIAN_FRONTEND=noninteractive apt install -y iproute2 cmake python python-pip build-essential gfortran wget curl
RUN pip install supervisor awscli
#RUN apt-get install openmpi-bin openmpi-common libopenmpi-dev iputils-ping -y
#ADD hordesat hordesat
ADD manyglucose-4.1-60 manyglucose-4.1-60

#RUN wget https://baldur.iti.kit.edu/hordesat/files/hordesat.zip
#RUN unzip hordesat.zip
RUN cd manyglucose-4.1-60/parallel && make
#ENV LD_LIBRARY_PATH=/usr/lib/openmpi/lib/:$LD_LIBRARY_PATH
#ADD test.cnf supervised-scripts/test.cnf
#ADD mpi-run.sh supervised-scripts/mpi-run.sh
ADD run.sh supervised-scripts/run.sh
#ADD make_combined_hostfile.py supervised-scripts/make_combined_hostfile.py
#RUN chmod 755 supervised-scripts/mpi-run.sh
RUN chmod 755 supervised-scripts/run.sh
#EXPOSE 22

#CMD hordesat/hordesat
#CMD supervised-scripts/mpi-run.sh
CMD supervised-scripts/run.sh
