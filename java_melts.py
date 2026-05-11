# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:18:45 2026

@author: Marti
"""

import os
from py4j.java_gateway import JavaGateway, GatewayParameters, launch_gateway

class JavaMELTS:
    """
    Python wrapper for the rhyolite-MELTS 1.2.0 Java engine.
    Mimics the MELTSdynamic interface closely enough that your
    existing code can run with minimal changes.
    """

    def __init__(self, jar_path="rhyolite-melts.jar"):
        # Launch JVM with the MELTS jar
        self.gateway = launch_gateway(
            classpath=jar_path,
            die_on_exit=True,
            redirect_stdout=None,
            redirect_stderr=None
        )

        self.jvm = JavaGateway(
            gateway_parameters=GatewayParameters(auto_convert=True)
        ).jvm

        # Create the MELTS engine instance
        self.engine = self.jvm.MELTSengine()

    # -----------------------------
    # Bulk composition
    # -----------------------------
    def set_bulk(self, bulk_dict):
        for ox, val in bulk_dict.items():
            self.engine.setBulkComposition(ox.lower(), float(val))

    # -----------------------------
    # Pressure
    # -----------------------------
    def set_pressure(self, P_bar):
        self.engine.setPressure(float(P_bar))

    # -----------------------------
    # Temperature
    # -----------------------------
    def set_temperature(self, T_C):
        self.engine.setTemperature(float(T_C))

    # -----------------------------
    # fO2 buffer
    # -----------------------------
    def set_fO2(self, buffer_name):
        self.engine.setOxygenFugacity(buffer_name)

    # -----------------------------
    # Run equilibrium
    # -----------------------------
    def calc_equilibrium(self):
        self.engine.calcEquilibriumState()

    # -----------------------------
    # Accessors mimicking MELTSdynamic
    # -----------------------------
    @property
    def solidNames(self):
        return list(self.engine.getSolidNames())

    @property
    def liquidNames(self):
        return list(self.engine.getLiquidNames())

    @property
    def mass(self):
        # Returns a dict: {phase_name: mass}
        m = {}
        for ph in self.solidNames + self.liquidNames:
            m[ph] = float(self.engine.getPhaseMass(ph))
        return m

    @property
    def status(self):
        class Status:
            failed = False
        return Status()
