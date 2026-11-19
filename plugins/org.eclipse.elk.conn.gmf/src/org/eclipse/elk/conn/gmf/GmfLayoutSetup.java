/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.conn.gmf;

import java.util.Collection;

import org.eclipse.core.runtime.IAdaptable;
import org.eclipse.elk.core.service.IDiagramLayoutConnector;
import org.eclipse.elk.core.service.ILayoutConfigurationStore;
import org.eclipse.elk.core.service.ILayoutSetup;
import org.eclipse.gmf.runtime.diagram.ui.editparts.IGraphicalEditPart;
import org.eclipse.gmf.runtime.diagram.ui.parts.DiagramEditor;

import com.google.inject.Binder;
import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Module;
import com.google.inject.util.Modules;

/**
 * Layout setup for the generic GMF connector.
 */
public class GmfLayoutSetup implements ILayoutSetup {

    @Override
    public boolean supports(final Object object) {
        if (object instanceof Collection) {
            Collection<?> collection = (Collection<?>) object;
            for (Object o : collection) {
                if (o instanceof IGraphicalEditPart) {
                    return true;
                }
            }
            return false;
        }
        return object instanceof DiagramEditor || object instanceof IGraphicalEditPart
                || object instanceof IAdaptable && ((IAdaptable) object).getAdapter(DiagramEditor.class) != null;
    }
    
    @Override
    public Injector createInjector(final Module defaultModule) {
        return Guice.createInjector(Modules.override(defaultModule).with(new GmfLayoutModule()));
    }
    
    /**
     * Guice module for the generic GMF connector.
     */
    public static class GmfLayoutModule implements Module {

        @Override
        public void configure(final Binder binder) {
            binder.bind(IDiagramLayoutConnector.class).to(GmfDiagramLayoutConnector.class);
            binder.bind(ILayoutConfigurationStore.Provider.class).to(GmfLayoutConfigurationStore.Provider.class);
        }
        
    }

}
