/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot;


import org.eclipse.elk.alg.graphviz.dot.serializer.GraphvizDotSerializer;
import org.eclipse.xtext.serializer.ISerializer;

/**
 * Use this class to register components to be used at runtime / without the Equinox extension registry.
 * 
 * @author msp
 */
public class GraphvizDotRuntimeModule extends AbstractGraphvizDotRuntimeModule {
    
    @Override
    public Class<? extends ISerializer> bindISerializer() {
        return GraphvizDotSerializer.class;
    }
    
}
