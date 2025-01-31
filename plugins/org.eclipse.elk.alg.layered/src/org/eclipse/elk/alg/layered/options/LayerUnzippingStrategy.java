/*******************************************************************************
 * Copyright (c) 2024 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 * 
 * SPDX-License-Identifier: EPL-2.0 
 *******************************************************************************/
package org.eclipse.elk.alg.layered.options;

/**
 * Strategies for unzipping layers i.e. splitting up nodes into multiple layers to
 * create more compact drawings.
 *
 */
public enum LayerUnzippingStrategy {
    
    NONE,
    /** Splits all layers with more than two nodes into several layers in an alternating pattern. */
    ALTERNATING;

}
