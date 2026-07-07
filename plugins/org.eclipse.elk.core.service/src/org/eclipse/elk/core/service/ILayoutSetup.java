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
package org.eclipse.elk.core.service;

import com.google.inject.Injector;
import com.google.inject.Module;

/**
 * Interface for Guice dependency injection setups used in the context of connecting ELK to
 * diagram viewers or editors. Such a setup must define at least a binding for {@link IDiagramLayoutConnector}.
 * A binding for {@link ILayoutConfigurationStore} can optionally be provided for enabling the Layout view.
 */
public interface ILayoutSetup {
    
    /**
     * Determine whether the layout connector provided by this layout setup is able to handle the given object.
     * 
     * @param object a workbench part or selected element, e.g. an edit part
     * @return true if this layout setup supports the object
     */
    boolean supports(Object object);
    
    /**
     * Create a Guice injector binding all implementations for this setup.
     * 
     * @param defaultModule a module with default bindings that should be considered, but may be overridden
     * @return a Guice injector
     */
    Injector createInjector(Module defaultModule);

}
