/*******************************************************************************
 * Copyright (c) 2018 TypeFox GmbH and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.graph.text.ui.layout;

import org.eclipse.elk.core.service.IDiagramLayoutConnector;
import org.eclipse.elk.core.service.ILayoutSetup;
import org.eclipse.elk.graph.text.ui.ElkGraphTextUiActivator;
import org.eclipse.xtext.ui.editor.XtextEditor;

import com.google.inject.Binder;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Module;
import com.google.inject.util.Modules;

/**
 * Setup for the {@link TextLayoutConnector}.
 */
public class TextLayoutSetup implements ILayoutSetup {

    @Override
    public boolean supports(Object object) {
        return object instanceof XtextEditor && ((XtextEditor) object).getLanguageName()
                .equals(ElkGraphTextUiActivator.ORG_ECLIPSE_ELK_GRAPH_TEXT_ELKGRAPH);
    }

    @Override
    public Injector createInjector(final Module defaultModule) {
        return Guice.createInjector(Modules.override(defaultModule).with(new TextLayoutModule()));
    }
    
    public static class TextLayoutModule implements Module {

        @Override
        public void configure(final Binder binder) {
            binder.bind(IDiagramLayoutConnector.class).to(TextLayoutConnector.class);
        }
        
    }

}
