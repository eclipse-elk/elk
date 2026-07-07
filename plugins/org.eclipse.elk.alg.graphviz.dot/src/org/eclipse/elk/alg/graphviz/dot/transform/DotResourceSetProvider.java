/*******************************************************************************
 * Copyright (c) 2011, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.alg.graphviz.dot.transform;

import org.eclipse.emf.ecore.resource.ResourceSet;
import org.eclipse.xtext.resource.XtextResourceSet;

import com.google.inject.Inject;
import com.google.inject.Provider;

/**
 * A graph format handler for Graphviz Dot. Instantiate through a proper injector:
 * 
 * <code>
 * new GraphvizDotStandaloneSetup()
 *     .createInjectorAndDoEMFRegistration()
 *     .getInstance(DotResourceSetProvider.class);
 * </code>
 *
 * @author cds
 */
public final class DotResourceSetProvider {
    
    @Inject private Provider<XtextResourceSet> resourceSetProvider;
    
    /**
     * Returns a new resource set for dot files.
     * 
     * @return resource set.
     */
    public ResourceSet createResourceSet() {
        return resourceSetProvider.get();
    }
    
}
