/*******************************************************************************
 * Copyright (c) 2016 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.service.internal;

import org.eclipse.elk.core.IGraphLayoutEngine;
import org.eclipse.elk.core.RecursiveGraphLayoutEngine;
import org.eclipse.elk.core.service.LayoutConnectorsService;

import com.google.inject.Binder;
import com.google.inject.Module;

/**
 * Default bindings that are included in all injectors created by the {@link LayoutConnectorsService}.
 */
public class DefaultModule implements Module {

    @Override
    public void configure(final Binder binder) {
        binder.bind(IGraphLayoutEngine.class).to(RecursiveGraphLayoutEngine.class);
    }

}
